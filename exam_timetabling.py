import csp
import time
import pandas as pd
import sys

class TimeTable(csp.CSP):

    def __init__(self):
        self.variables = []
        self.domains = dict()
        self.neighbors = dict()
        self.weight = {}    #for dom/wdeg

        self.countOfVarChecks= 0
        self.days = 21
        self.timeSlots = 3
        self.df = pd.read_csv("Στοιχεία Μαθημάτων.csv")

        semesterColumn = self.df["Εξάμηνο"]
        subjectColumn = self.df["Μάθημα"]
        labsColumn = self.df["Εργαστήριο (TRUE/FALSE)"]
        difficultColumn = self.df["Δύσκολο (TRUE/FALSE)"]
        professorColumn = self.df["Καθηγητής"]

        #variables -> the titles of the subjects
        #     |Day 1 | Day 2 | .... | Day 21
        #09-12|  X1      X3             Xm-1
        #12-15|  X2      --             Xm
        #15-18|  --      X4             --

        labNum = 0          #we create pairs of theory-lab subjects using this counter
        self.varsInfo = {}  #store usefull info about each subject

        for subject, labSubject, semester, difficult, prof in zip(subjectColumn.values, labsColumn.values, 
        semesterColumn.values, difficultColumn.values, professorColumn.values):
            #if lab subject exists we add the labnum and also we create a new subject for the lab
            if labSubject:
                self.varsInfo[subject] = (subject, semester, labNum, difficult, prof)   #we store the info
                self.variables.append(subject)  #we update the var list

                labsubject = subject + " Lab"
                self.varsInfo[labsubject] = (labsubject, semester, labNum, difficult, prof)
                self.variables.append(labsubject) #we add the labs as subjects also

                labNum += 1

            #if the subject doesnt have a lab we put None as labnum
            else:
                self.variables.append(subject) #we update the var list
                self.varsInfo[subject] = (subject, semester, None, difficult, prof)

        #domains -> the place in the table
        #     |Day 1 | Day 2 | .... | Day 21
        #09-12|   1      4             
        #12-15|   2      5
        #15-18|   3      6             3*21=63
        for neighboor in self.variables:
            var = self.varsInfo[neighboor]
            currDomain = []
            for i in range(1, self.days*self.timeSlots + 1):
                if "Lab" in var[0] and i%3 == 1: continue                              #labs cant be in the first slot of the day(=> i%3==1 )
                elif var[2] is not None and "Lab" not in var[0] and i%3 == 0: continue #theory with lab cant be in the last slot of the day(=> i%3==0)
                else: currDomain.append(i)
            self.domains[neighboor] = currDomain

        #neighboors -> each subject has all the other subjects as neighboors (expect itself)
        for neighboor in self.variables:
            self.neighbors[neighboor] = [x for x in self.variables if x!=neighboor]
        
        #for dom/wdeg
        makeWeights(self)
        
        #we create the csp
        csp.CSP.__init__(self, self.variables, self.domains, self.neighbors, self.var_constraints)

    # -------------------------------------- Constraints ---------------------------------------------------------#
    def semesterConstraint(self, A, a, B, b):
        #we check that the different subjects(labs excluded) of the day are from different semesters
        #we check to see if the semester is the same
        if A[1] == B[1]:
            #we check if the relation between A and B is not a theory-lab
            if not(A[2] is not None and A[2] == B[2]):
                #if the min is on the fist time period 9-12 the other must be atleast 3 slots away(until the next day)
                if min(a,b)%3 == 1 and abs(a-b) < 3: return False
                #or if the mins is on the second time period 12-15 the must be atleast 2 slots away
                elif min(a,b)%3 == 2 and abs(a-b) < 2: return False
        return True

    def theoryLabConstraint(self, A, a, B, b):
        #if A and B are a pair(subject-lab) we check if they are next to each other
        if A[2] is not None and A[2] == B[2]:
            if "Lab" in A[0]:
                # the slot of the lab should be right after the slot of the theory
                if (b > a) or (abs(a-b) != 1): return False
            else:
                if (a > b) or (abs(a-b) != 1): return False  
        return True

    def difficultyConstraint(self, A, a, B, b):
        #we check if difficult subjects are too close (labs excluded)
        if A[3] and B[3]:
            #every 2 difficult subjects must have 1 day difference
            if not(A[2] is not None and A[2] == B[2]):
                #so if the min is on the first time slot the other must be at least 6 slots away 
                #(2 slot of the same day and other 3 for the empty day)
                if ((min(a,b)%3 == 1 and abs(a-b) < 6)

                #if the min is on the second time slot the other must be at least 5 slots away 
                #(1 slot of the same day and other 3 for the empty day)
                or (min(a,b)%3 == 2 and abs(a-b) < 5)

                #if the min is on the third time slot the other must be at least 4 slots away 
                #3(for the empty day)
                or (min(a,b)%3 == 0 and abs(a-b) < 4)): return False
        return True

    def professorConstraint(self, A, a, B, b):
        #we check if the same professor has multiple subjects in the same day (labs excluded)
        if A[4] == B[4]:
            if not(A[2] is not None and A[2] == B[2]):
                if (min(a,b)%3 == 1 and abs(a-b) < 3
                or min(a,b)%3 == 2 and abs(a-b) < 2): return False
        return True
    # -----------------------------------------------------------------------------------------------------------#

    def var_constraints(self, A, a, B, b):
        #we count the number of time we need to check the constraints
        self.countOfVarChecks += 1

        # we find the subjects info
        A = self.varsInfo[A]
        B = self.varsInfo[B]

        #only check that on each time-slot we have one subject
        if a == b: return False

        if(not self.semesterConstraint(A,a,B,b)): return False

        if(not self.theoryLabConstraint(A,a,B,b)): return False

        if(not self.difficultyConstraint(A,a,B,b)): return False

        if(not self.professorConstraint(A,a,B,b)): return False

        #if we pass all the constraints we are good to go :)
        return True

    #for the index :
    #0-> prints the names, 1-> prints the semesters 2-> peints the lab-theory pairs
    #3-> the difficulty 4-> the professor
    def printTimeTalbe(self,solution,index):

        #You can adjust the {:80} for your display preferences
        #-----------------------------------------------------

        #we must sort first the solution by value to print out correctly
        solution = dict(sorted(solution.items(), key=lambda item: item[1]))

        maxSubsPerDay = 0
        day = 1
        counter = 1

        print('{:6s}'.format("\nDay" + str(day)), end = " ")
        for subjectInfo, timeslot in solution.items():
            curr = self.varsInfo[subjectInfo]
            if counter == timeslot:
                print('{:2d} {:80}'.format(timeslot,str(curr[index])), end = " ")
                counter += 1
                maxSubsPerDay += 1
                if maxSubsPerDay >= 3 :
                    day += 1
                    maxSubsPerDay = 0
                    print('{:6s}'.format("\nDay" + str(day)), end = " ")
            else:
                for emptySlots in range(timeslot-counter):
                    print('{:2d} {:80}'.format(counter, "----------"), end = " ")
                    counter += 1
                    maxSubsPerDay += 1
                    if maxSubsPerDay >= 3 :
                        day += 1
                        maxSubsPerDay = 0
                        print('{:6s}'.format("\nDay" + str(day)), end = " ")
                print('{:2d} {:80}'.format(timeslot,str(curr[index])), end = " ")
                maxSubsPerDay += 1
                if maxSubsPerDay >= 3 :
                    day += 1
                    maxSubsPerDay = 0
                    print('{:6s}'.format("\nDay" + str(day)), end = " ")  
                counter = timeslot + 1
        print('\n')

#initialize weight for dom/wdeg
def makeWeights(csp):
    for var in csp.variables:
        for n in csp.neighbors[var]:
            csp.weight[(var,n)] = 1

def AC3(csp, queue=None, removals=None, arc_heuristic=csp.dom_j_up):
    """[Figure 6.3]"""
    if queue is None:
        queue = {(Xi, Xk) for Xi in csp.variables for Xk in csp.neighbors[Xi]}
    csp.support_pruning()
    queue = arc_heuristic(csp, queue)
    checks = 0
    while queue:
        (Xi, Xj) = queue.pop()
        revised, checks = revise(csp, Xi, Xj, removals, checks)
        if revised:
            if not csp.curr_domains[Xi]:
                return False, checks  # CSP is inconsistent
            for Xk in csp.neighbors[Xi]:
                if Xk != Xj:
                    queue.add((Xk, Xi))
    return True, checks  # CSP is satisfiable

def revise(csp, Xi, Xj, removals, checks=0):
    """Return true if we remove a value."""
    revised = False
    for x in csp.curr_domains[Xi][:]:
        # If Xi=x conflicts with Xj=y for every possible y, eliminate Xi=x
        # if all(not csp.constraints(Xi, x, Xj, y) for y in csp.curr_domains[Xj]):
        conflict = True
        for y in csp.curr_domains[Xj]:
            if csp.constraints(Xi, x, Xj, y):
                conflict = False
            checks += 1
            if not conflict:
                break
        if conflict:
            csp.prune(Xi, x, removals)
            revised = True
    if not csp.curr_domains[Xi]:
        csp.weight[(Xi,Xj)] += 1
        csp.weight[(Xj,Xi)] += 1
    return revised, checks

def mac(csp, var, value, assignment, removals, constraint_propagation=AC3):
    """Maintain arc consistency."""
    return constraint_propagation(csp, {(X, var) for X in csp.neighbors[var]}, removals)

def forward_checking(csp, var, value, assignment, removals):
    """Prune neighbor values inconsistent with var=value."""
    csp.support_pruning()
    for B in csp.neighbors[var]:
        if B not in assignment:
            for b in csp.curr_domains[B][:]:
                if not csp.constraints(var, value, B, b):
                    csp.prune(B, b, removals)
            if not csp.curr_domains[B]:
                csp.weight[(var,B)] += 1
                csp.weight[(B,var)] += 1
                return False
    return True

def dom_wdeg(assignment, csp):
    minVal = float('inf')
    bestVar = None

    for var in csp.variables:
        if var in assignment:
            continue

        dom = len(csp.choices(var))
        wdeg = 1
        for n in csp.neighbors[var]:
            wdeg += csp.weight[(var,n)]
        
        ratio = (float)(dom/wdeg)
        if ratio < minVal:
            minVal = ratio
            bestVar = var
    return bestVar

#function just to get the terminal arguments
def getArguments():
    variableOrdering = None
    if(sys.argv[1] == "first_unassigned_variable"): 
        print("Variable ordering: first_unassigned_variable")
        variableOrdering = csp.first_unassigned_variable
    elif(sys.argv[1] == "mrv"):
        print("Variable ordering: mrv")
        variableOrdering = csp.mrv
    elif(sys.argv[1] == "dom_wdeg"):
        print("Variable ordering: dom_wdeg")
        variableOrdering = dom_wdeg
    else:
        print("Wrong Variable ordering")
    
    valueOrdering = None
    if(sys.argv[2] == "unordered_domain_values"): 
        print("Value ordering: unordered_domain_values")
        valueOrdering = csp.unordered_domain_values
    elif(sys.argv[2] == "lcv"):
        print("Value ordering: lcv")
        valueOrdering = csp.lcv
    else:
        print("Wrong Value ordering")

    inference = None
    if(sys.argv[3] == "no_inference"): 
        print("Inference: no_inference")
        inference = csp.no_inference
    elif(sys.argv[3] == "forward_checking"):
        print("Inference: forward_checking")
        if sys.argv[1] == "dom_wdeg":
            inference = forward_checking
        else:
            inference = csp.forward_checking
    elif(sys.argv[3] == "mac"):
        print("Inference: mac")
        inference = csp.mac
        if sys.argv[1] == "dom_wdeg":
            inference = mac
        else:
            inference = csp.mac
    else:
        print("Wrong Inference")

    return(variableOrdering,valueOrdering,inference)

if __name__ == '__main__':

    problem = TimeTable()
    solution = None
    begin = None
    end = None 

    #without args we run the min conflicts
    if len(sys.argv) == 1:
        print("Inference: min_conflicts")
        begin = time.time()
        solution = csp.min_conflicts(problem)
        end = time.time()
    
    #with args we run the backtracking_search
    else:
        args = getArguments()
        #we check if the input was correct
        for i in range(3):
            if args[i] is None: exit(1)

        begin = time.time()
        solution = csp.backtracking_search(problem, args[0], args[1] , args[2])
        end = time.time()

    TimeTable.printTimeTalbe(problem,solution,0)
    print("The total time is: \n" + str(end - begin))
    print("The total assigments are: " + str(problem.nassigns))
    print("The total constraints checks are:",problem.countOfVarChecks)

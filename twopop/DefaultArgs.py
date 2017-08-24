import sys

def DefaultArgs(givenArgs, defArgs):
    nArgs = len(defArgs)
    out = []
    if(nArgs != len(givenArgs)):
        print "inside DefaultArgs : assigning default args" 
        nGivenArgs = len(givenArgs)
        argCount = 1
        for ll, lArg in enumerate(defArgs):
            if(argCount > nGivenArgs or givenArgs[ll] == '[]'):
                out.append(lArg)
            else:
                out.append(givenArgs[ll])
            argCount = argCount + 1
    else:
        for ll, lArg in enumerate(givenArgs):
            if(lArg != '[]'):
                out.append(lArg)
            else :
                out.append(defArgs[ll])
    return out

    

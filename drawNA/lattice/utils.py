from bisect import bisect_left, bisect_right

def find_crossover_locations(
    max_size: int = 30,
    bp_per_turn: float = 10.45,
    kind: str = "half",
    origin: int = 0
) -> dict:
    """
    Function which returns an array containing the best locations 
    for half or whole turns of the DNA strand where the strand will
    stay in plane and have the least strain/stress build-up.
    
    This is computed given the number of basepairs per turn
    & locations are given in no. of base pairs

    Arguments:
    bp_per_turn - (default: 10.45)
    kind - either half turns or whole turns
    max_size - max no. of turns to base pairs to the list
    origin - what number crossovers start at, either 0 (python indexing) or 1 
    """

    no_of_bp = 0
    if origin == 0: # i.e. Python indexing
        crossover_location = [0]
        shift = 1
    else:
        crossover_location = [1]
        shift = 0


    # half turn
    if kind == "half":
        i = 1  # every odd number
        while no_of_bp < max_size:
            # minus 1 at the end because Python indexing starts at 0
            no_of_bp = int(round(bp_per_turn * 0.5 * i, 0) - shift)
            crossover_location.append(no_of_bp)
            i += 2

    # whole turns
    elif kind == "whole":
        i = 2  # every even number
        while no_of_bp < max_size:
            # minus 1 at the end because Python indexing starts at 0
            no_of_bp = int(round(bp_per_turn * 0.5 * i, 0) - shift)
            crossover_location.append(no_of_bp)
            i += 2
    return crossover_location

def find_closest(myList, myNumber):
    """
    Credit: https://stackoverflow.com/questions/12141150
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after + 1
    else:
        return before + 1
from drawNA.oxdna import Strand, Nucleotide, System
from drawNA.oxdna.strand import generate_helix, POS_BACK
import random
import numpy as np

FENE_LENGTH = 0.76

#Generation of initial strand configuration
def generate_system(box, length=16, n_strands=10, stapled=5, percentage=30):
    # initially create a strand of length n that is double stranded
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True, start_position=np.array([0., 0., 0.]))
    strands.append(strand.copy())
    doubles.append(double.copy())

    if n_strands > 1:
        for i in range(n_strands-1):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1

            # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand, double = generate_helix(
                n=length,
                start_position=start,
                direction=direction,
                a1=a1,
                double=True,
            )
            strands.append(strand)
            doubles.append(double)


        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # add an additional strand with only 46 bases which is the rest of the M13 strand to be simulated
        strand, double = generate_helix(
            n=50,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)


         #   # using the two previously created strands to create a new strand that is added to the system
         #   nucleotides = []
         #   for strand in strands:
         #       nucleotides += strand.nucleotides

         #   strand = Strand(nucleotides=nucleotides)

         #   # create system and add the final completed strand
         #   main_system = System(box)
         #   main_system.add_strand(strand)

         #   actual_doubles = []
         #   for strand in doubles:

         #       nucleotides = strand.nucleotides[:stapled] #append only a certain number of doubles (actual doubles)
         #       actual_doubles.append(Strand(nucleotides=nucleotides))

         #   main_system.add_strands(actual_doubles)


    # using the two previously created strands to create a new strand that is added to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides

    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    main_system = System(box)
    main_system.add_strand(strand)


    # pick a random number between 0 and stapled (e.g. 10 bases)
    # keep picking random number until certain percentage is reached
    actual_doubles = []
    for strand in doubles[:-1]:
        block_number = length/stapled
        ds_block_number = block_number * (percentage*0.01)
        new_staples = random.sample(range(int(block_number)), int(ds_block_number))
        for block in new_staples:
            base_start = block * stapled
            base_end = base_start + stapled
            nucleotides = strand.nucleotides[base_start:base_end] # append only a certain number of doubles (actual doubles)
            actual_doubles.append(Strand(nucleotides=nucleotides))

    for strand in doubles[-1:]:
        block_number = (50)/stapled
        ds_block_number = block_number * ((percentage)*0.01)
        new_staples = random.sample(range(int(block_number)), int(ds_block_number))
        for block in new_staples:
            base_start = block * stapled
            base_end = base_start + stapled
            nucleotides = strand.nucleotides[base_start:base_end] # append only a certain number of doubles (actual doubles)
            actual_doubles.append(Strand(nucleotides=nucleotides))

    main_system.add_strands(actual_doubles)

    return main_system


#Generation of initial strand configuration
def generate_system_all_random(box, length=16, n_strands=10, stapled=5, percentage=30):
    """
    generate a folded strand with randomly distributed staples along the WHOLE strand
    """
    # initially create a strand of length n that is double stranded
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True, start_position=np.array([0., 0., 0.]))
    strands.append(strand.copy())
    doubles.append(double.copy())

    if n_strands > 1:
        for i in range(n_strands-1):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1

            # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand, double = generate_helix(
                n=length,
                start_position=start,
                direction=direction,
                a1=a1,
                double=True,
            )
            strands.append(strand)
            doubles.append(double)

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # add an additional strand with only 46 bases which is the rest of the M13 strand to be simulated
        strand, double = generate_helix(
            n=64,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)


    # using the two previously created strands to create a new strand that is added to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides

    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    main_system = System(box)
    main_system.add_strand(strand)

    actual_doubles = []

    block_number = (8070)/stapled
    ds_block_number = block_number * (percentage*0.01)
    new_staples = random.sample(range(int(block_number)), int(ds_block_number))

    # there are 40 blocks per row but 5 in the last entry of doubles
    for starting_point in new_staples:
        print(starting_point)
        hairpin = int(starting_point // (length/stapled))
        print(hairpin)

        strand = doubles[hairpin]
        print(strand)

        base_start_beginning = int((starting_point - (hairpin*(length/10)))*10)

        base_end_beginning = base_start_beginning + stapled
        if base_end_beginning > len(strand):
            base_end_beginning = len(strand)

        base_start = len(strand) - base_end_beginning
        base_end = len(strand) - base_start_beginning


        # append only a certain number of doubles (actual doubles)
        nucleotides = strand.nucleotides[base_start:base_end]
        actual_doubles.append(Strand(nucleotides=nucleotides))

    main_system.add_strands(actual_doubles)

    return main_system

def generate_loop_random(box, length=16, n_strands=10, stapled=5, percentage=30):
    """
    generate a folded strand with randomly distributed staples along the WHOLE strand
    """
    # initially create a strand of length n that is double stranded
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True, start_position=np.array([0., 0., 0.]))
    strands.append(strand.copy())
    doubles.append(double.copy())

    if n_strands > 1:
        #for i in range(n_strands-1):
        for i in range(int((n_strands)/2)-1):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1
            # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand, double = generate_helix(
                n=length,
                start_position=start,
                direction=direction,
                a1=a1,
                double=True,
            )
            strands.append(strand)
            doubles.append(double)

        #add a strand that goes in the same direction than the previous one to initiate loop
        last_nuc = strands[-1].nucleotides[-1]
        direction = last_nuc._a3
        a1 = last_nuc._a1
        start = last_nuc.pos_back + [-4e-01, -5e-01, -1e-01]#+ (FENE_LENGTH - POS_BACK) * a1

        strand, double = generate_helix(
            n=432,#164,#188
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)

        for i in range(  int(((n_strands)/2)), int(n_strands-1)  ):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1
            # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand, double = generate_helix(
                n=length,
                start_position=start,
                direction=direction,
                a1=a1,
                double=True,
            )
            strands.append(strand)
            doubles.append(double)

        #add a strand that goes in the same direction than the previous one to initiate loop
        last_nuc = strands[-1].nucleotides[-1]
        direction = last_nuc._a3
        a1 = last_nuc._a1
        start = last_nuc.pos_back + [4e-01, 5e-01, -1e-01] #+ (FENE_LENGTH - POS_BACK) * a1
        print(start)

        strand, double = generate_helix(
            n=32,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)
        #print(strand)




#        last_nuc = strands[-1].nucleotides[-1]
#        direction = -last_nuc._a3
#        a1 = -last_nuc._a1

#        # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
#        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1
#
#        # add an additional strand with only 46 bases which is the rest of the M13 strand to be simulated
#        strand, double = generate_helix(
#            n=64,
#            start_position=start,
#            direction=direction,
#            a1=a1,
#            double=True,
#        )
#        strands.append(strand)
#        doubles.append(double)


    # using the two previously created strands to create a new strand that is added to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides


    strand = Strand(nucleotides=nucleotides)

    #print(strand.nucleotides[0], strand.nucleotides[-1], strand.nucleotides[-1].index, strand.nucleotides[0].index )
    strand.circular = True
    #print(strand.nucleotides[0], strand.nucleotides[-1], strand.nucleotides[-1].index, strand.nucleotides[0].index )

    # create system and add the final completed strand
    main_system = System(box)
    main_system.add_strand(strand)

    print(main_system, strand)

    actual_doubles = []

    block_number = (8070)/stapled 
    #block_number = (6530)/stapled
    ds_block_number = block_number * (percentage*0.01)
    new_staples = random.sample(range(int(block_number)), int(ds_block_number))

    # there are 40 blocks per row but 5 in the last entry of doubles
    for starting_point in new_staples:
        print(starting_point)
        hairpin = int(starting_point // (length/stapled))
        print(hairpin)

        strand = doubles[hairpin]
        print(strand)

        base_start_beginning = int((starting_point - (hairpin*(length/10)))*10)

        base_end_beginning = base_start_beginning + stapled
        if base_end_beginning > len(strand):
            base_end_beginning = len(strand)

        base_start = len(strand) - base_end_beginning
        base_end = len(strand) - base_start_beginning


        # append only a certain number of doubles (actual doubles)
        nucleotides = strand.nucleotides[base_start:base_end]
        actual_doubles.append(Strand(nucleotides=nucleotides))

    main_system.add_strands(actual_doubles)

    

    return main_system



#Generation of initial strand configuration
def generate_system_single(box, length=16, n_strands=10, stapled=5, percentage=30):
    # initially create a strand of length n that is double stranded
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True, start_position=np.array([0., 0., 0.]))
    strands.append(strand.copy())
    doubles.append(double.copy())

    if n_strands > 1:
        for i in range(n_strands-1):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1

            # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand, double = generate_helix(
                n=length,
                start_position=start,
                direction=direction,
                a1=a1,
                double=True,
            )
            strands.append(strand)
            doubles.append(double)

###### this is needed if we do 8050 bases!!!
#        last_nuc = strands[-1].nucleotides[-1]
#        direction = -last_nuc._a3
#        a1 = -last_nuc._a1

#        # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
#        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

#        # add an additional strand with only 46 bases which is the rest of the M13 strand to be simulated
#        strand, double = generate_helix(
#            n=50,
#            start_position=start,
#            direction=direction,
#            a1=a1,
#            double=True,
#        )
#        strands.append(strand)
#        doubles.append(double)



         #   # using the two previously created strands to create a new strand that is added to the system
         #   nucleotides = []
         #   for strand in strands:
         #       nucleotides += strand.nucleotides

         #   strand = Strand(nucleotides=nucleotides)

         #   # create system and add the final completed strand
         #   main_system = System(box)
         #   main_system.add_strand(strand)

         #   actual_doubles = []
         #   for strand in doubles:

         #       nucleotides = strand.nucleotides[:stapled] #append only a certain number of doubles (actual doubles)
         #       actual_doubles.append(Strand(nucleotides=nucleotides))

         #   main_system.add_strands(actual_doubles)


    # using the two previously created strands to create a new strand that is added to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides

    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    main_system = System(box)
    main_system.add_strand(strand)


    # pick a random number between 0 and stapled (e.g. 10 bases)
    # keep picking random number until certain percentage in reached

  #  actual_doubles = []
  #  for strand in doubles:
  #      block_number = length/stapled
  #      base = 0
  #      for block in range block_number:
  #          nucleotides = strand.nucleotides[base:stapled] #append only a certain number of doubles (actual doubles)
  #          actual_doubles.append(Strand(nucleotides=nucleotides))
  #          base += stapled


    actual_doubles = []
  #  for strand in doubles[:-1]:
    for strand in doubles:
        block_number = length/stapled
        ds_block_number = block_number * (percentage*0.01)
        #print(block_number, ds_block_number)

        new_staples = random.sample(range(int(block_number)), int(ds_block_number))
       #print(new_staples)

        for block in new_staples:
            base_start = block * stapled
            base_end = base_start + stapled
            nucleotides = strand.nucleotides[base_start:base_end] # append only a certain number of doubles (actual doubles)
            actual_doubles.append(Strand(nucleotides=nucleotides))
            #print(actual_doubles)

###### this is needed if we do 8050 bases!!!  
#    for strand in doubles[-1:]:
#        block_number = (50)/stapled
#        ds_block_number = block_number * ((percentage/2)*0.01)
#       # print(block_number, ds_block_number)#
#
#        new_staples = random.sample(range(int(block_number)), int(ds_block_number))
#       #print(new_staples)#
#
#        for block in new_staples:
#            base_start = block * stapled
#            base_end = base_start + stapled
#            nucleotides = strand.nucleotides[base_start:base_end] # append only a certain number of doubles (actual doubles)
#            actual_doubles.append(Strand(nucleotides=nucleotides))
#            #print(actual_doubles)



    main_system.add_strands(actual_doubles)

    return main_system


def generate_exact_system(box, length=16, n_strands=10, stapled=5, percentage=30, binary_file='20_string.txt'):
    """
    generate a system with a given binary input file (0=ss, 1=ds)
    """

    # initially create a strand of length n that is double stranded
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True, start_position=np.array([0., 0., 0.]))
    strands.append(strand.copy())
    doubles.append(double.copy())

    if n_strands > 1:
        for i in range(n_strands-1):

            last_nuc = strands[-1].nucleotides[-1]
            direction = -last_nuc._a3
            a1 = -last_nuc._a1

            # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
            start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

            # generate strand above that's going in opposite direction
            strand, double = generate_helix(
                n=length,
                start_position=start,
                direction=direction,
                a1=a1,
                double=True,
            )
            strands.append(strand)
            doubles.append(double)

        
        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # add an additional strand with only 46 bases which is the rest of the M13 strand to be simulated
        strand, double = generate_helix(
            n=64,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)

    # using the two previously created strands to create a new strand that is added to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides

    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    main_system = System(box)
    main_system.add_strand(strand)


    # add only the doubles
    with open(binary_file, 'r') as binary_system:
        DNA_system = binary_system.read()

    #create a list new_staples = [] that contains all indices for each first '1' 
    new_staples = []
    new_staples_end = []

    #DNA_system.index('01')
    import re
    new_staples = [m.start() for m in re.finditer('01', DNA_system)]
    new_staples = [x+1 for x in new_staples]
    #new_staples.insert(0, 0)

    new_staples_end = [m.start() for m in re.finditer('10', DNA_system)]
    new_staples_end = [x+1 for x in new_staples_end]

    new_staples_end.append(8064)
    #print(len(new_staples), len(new_staples_end))

    actual_doubles = []

    index = 0
    for starting_point in new_staples:

        hairpin = int(starting_point // length)
        strand = doubles[hairpin]
        

        base_start_beginning = int((starting_point - (hairpin*length)))
        base_end_beginning = base_start_beginning + stapled

        if base_end_beginning > len(strand):
            base_end_beginning = len(strand)

        base_start = len(strand) - base_end_beginning
        base_end = len(strand) - base_start_beginning
       
        nucleotides = strand.nucleotides[base_start:base_end]

        actual_doubles.append(Strand(nucleotides=nucleotides))
        print('first staples:')
        print(starting_point, hairpin, base_start, base_end)
        #new_base_start = len(strand) - (int((new_staples_end[index] - (hairpin*length))))

     #   if base_start != new_base_start and starting_point == new_staples[-1]:
     #       base_start = new_base_start
     #   else:
        new_base_start = len(strand) - (int((new_staples_end[index] - (hairpin*length))))

        while base_start != new_base_start:

#            hairpin = int(new_staples_end[index] // length)

#            base_end = base_start
#
#            if base_end == len(strand):
#                base_end == 0 

#            strand = doubles[hairpin]

         #   base_start = new_base_start 

#            if new_base_start < 0:
#                base_start = 0
#                nucleotides = strand.nucleotides[base_start:base_end]
#                actual_doubles.append(Strand(nucleotides=nucleotides))
#                print('additional staples:')
#                print(starting_point, hairpin, base_start, base_end)

#                base_start = len(strand) + new_base_start # plus because - and - is plus
#                hairpin = hairpin + 1
#                strand = doubles[hairpin]

#                nucleotides = strand.nucleotides[base_start:base_end]
#                actual_doubles.append(Strand(nucleotides=nucleotides))
#                print('additional staples:')
#                print(starting_point, hairpin, base_start, base_end)

#            else:
            hairpin = int(new_staples_end[index] // length)

            base_end = base_start

            if base_end == len(strand):
                base_end == 0 

            strand = doubles[hairpin]


        #    elif new_base_start == 0:
        #        hairpin = hairpin - 1

        #    if base_start == length:
        #        base_start = 0


        #    base_start = new_base_start 

            if new_base_start < 0:
                base_start = 0
                hairpin = hairpin - 1
                strand = doubles[hairpin]

                nucleotides = strand.nucleotides[base_start:base_end]
                actual_doubles.append(Strand(nucleotides=nucleotides))
                print('additional staples:')
                print(starting_point, hairpin, base_start, base_end)

                hairpin = hairpin + 1
                strand = doubles[hairpin]
                new_base_start = len(strand) + new_base_start 
                base_end = len(strand)

            elif new_base_start == 0 and hairpin != 20 :
                hairpin = hairpin - 1
                strand = doubles[hairpin]

            base_start = new_base_start 

            nucleotides = strand.nucleotides[base_start:base_end]
            actual_doubles.append(Strand(nucleotides=nucleotides))
            print('additional staples:')
            print(starting_point, hairpin, base_start, base_end)

        index += 1

    print(actual_doubles)
    main_system.add_strands(actual_doubles)

    return main_system

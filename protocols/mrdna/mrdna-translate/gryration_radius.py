#!/usr/bin/python

import math 
import sys
import os
import numpy as np
from datetime import datetime

def calculate_Rg(filename):
    '''
    Calculates the Radius of Gyration (Rg) of DNA (.pdb format) 
    Returns the Rg integer value in Angstrom.
    '''
    coord = []#list()
    mass = []#list()
    rg_all = []#list()
    Structure = open(filename, 'r')
    frame_number = 0

#    for frame in Structure
#        #for line in Structure:
#        if frame.startswith('END'):
#            frame_number += 1

    for line in Structure:
        #if not line.startswith('ATOM'):

        if line.startswith('CRYST1'):
            continue
        elif line.startswith('ATOM'):
            try:
                line = line.split()
                x = float(line[6])
                y = float(line[7])
                z = float(line[8])
                coord.append([x, y, z])
                if line[2] == 'DNA':
                    mass.append(650.0) #average weight of one bp is 650 dalton - 2bp are represented by one bead in this model
                elif line[2] == 'NAS':
                    mass.append(650.0) 
                elif line[2] == 'O':
                    mass.append(650.0) 
            except:
                pass
        elif line.startswith('END'):
            xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
            tmass = sum(mass)
            rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
            mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
            rg = math.sqrt(rr / tmass-mm)
            rg_all.append(round(rg, 3))
           # print(rg_all)
        
   # xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
   # tmass = sum(mass)
   # rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
   # mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
   # rg = math.sqrt(rr / tmass-mm)
    #print(len(rg_all))

    rg = sum(rg_all[-200::])/len(rg_all[-200::]) #[round(rg_all, 3)]
    max_rg = [max(rg_all[-200::])]
    min_rg = [min(rg_all[-200::])]

    print(filename)
    print(rg_all, rg, max_rg, min_rg)

    return(rg_all, rg, max_rg, min_rg)




def plot_Rg(rg, max_rg, min_rg, percentage_ds_exact, percentage_ds_random):
    import plotly.graph_objects as go

    above_avg = max_rg - rg
    below_avg = rg - min_rg

    a_avg_random = [above_avg[0],above_avg[1],above_avg[2], above_avg[3], above_avg[4],  above_avg[6],above_avg[7], above_avg[9], above_avg[10], above_avg[11], above_avg[12], above_avg[13], above_avg[14], above_avg[15],above_avg[16], above_avg[17], above_avg[19], above_avg[20], above_avg[21], above_avg[22], above_avg[23], above_avg[25], above_avg[26], above_avg[27] ]
    b_avg_random = [below_avg[0],above_avg[1],above_avg[2], above_avg[3], below_avg[4], below_avg[6], below_avg[7], below_avg[9], below_avg[10], above_avg[11], above_avg[12], above_avg[13], above_avg[14], above_avg[15], above_avg[16], above_avg[17], above_avg[19], above_avg[20], above_avg[21], above_avg[22], above_avg[23], above_avg[25], above_avg[26], above_avg[27] ]

    a_avg = [ above_avg[8], above_avg[12], above_avg[18], above_avg[24] ]
    b_avg = [below_avg[8], below_avg[12], below_avg[18], below_avg[24] ]


    fig = go.Figure(data=go.Scatter(
            x=percentage_ds_random,
            y=[rg[0], rg[1], rg[2], rg[3], rg[4], rg[6], rg[7], rg[9], rg[10], rg[11], rg[13], rg[14], rg[15], rg[16], rg[17],  rg[19], rg[20], rg[21], rg[22], rg[23],rg[25], rg[26], rg[27]],
            name="random", opacity=0.6, mode='markers',
            error_y=dict(
                type='data',
                symmetric=False,
                array=a_avg_random,
                arrayminus=b_avg_random)
            ))

    fig.add_scatter(x=percentage_ds_exact,y=[ rg[5], rg[9], rg[14], rg[19]], mode="markers", name='exact', opacity=0.6,
        error_y=dict(
            type='data',
            symmetric=False,
            array=a_avg,#[0.1, 0.2, 0.1, 0.1]#,
            arrayminus=b_avg)
        )

    fig.update_layout(
        title="Gyration radii over percentage of hybridised sections in DNA strands",
        xaxis_title="Hybridised DNA (%)",
        yaxis_title="Gyration radius (A)",
        legend_title="Hybridised section distribution"
       # font=dict(
       #     family="Courier New, monospace",
        #    size=18,
        #    color="RebeccaPurple"
        #)
    )

    fig.write_html("gyration_radii_both6.html")
    fig.write_image("gyration_radii_both6.png")

    #fig.show()

def main(args):
    pdb_folder = args['folder']

    all_rg = np.array([])
    all_max_rg = np.array([])
    all_min_rg = np.array([])
    all_rg_all = np.array([])

    print('>>>>>starting collecting data from pdb files!')

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

    # for each pdb file calculate the gyration radius (including variance over frames) 
    # and append to array containing all gyration values for all pdb files in the folder

    for pdbfile in sorted(os.listdir(pdb_folder)):
        print(sorted(os.listdir(pdb_folder)))
        pdbfile_path = pdb_folder + '/' + pdbfile
        rg_all, rg, max_rg, min_rg = calculate_Rg(pdbfile_path)

        all_rg_all = np.append(all_rg_all, rg_all)
        all_rg = np.append(all_rg, rg)
        all_max_rg = np.append(all_max_rg, max_rg)
        all_min_rg = np.append(all_min_rg, min_rg)
    
    print(all_rg_all, all_rg, all_max_rg, all_min_rg)
    print('>>>>>finished collecting data from pdb files!')

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

    percentage_ds_exact = [20, 40, 60, 80]
    percentage_ds_random = [0, 0, 0, 0, 100, 100, 100, 100, 20, 20, 20, 20, 40, 40, 40, 40, 40, 50, 60, 60,60, 60, 60, 70, 80, 80, 80, 90]

    plot_Rg(all_rg, all_max_rg, all_min_rg, percentage_ds_exact, percentage_ds_random)

    np.savetxt("all_rg.csv", all_rg, delimiter=",")
    np.savetxt("all_max_rg.csv", all_max_rg, delimiter=",")
    np.savetxt("all_min_rg.csv", all_min_rg, delimiter=",")
    np.savetxt("percentage_ds_exact.csv", percentage_ds_exact, delimiter=",")
    np.savetxt("percentage_ds_random.csv", percentage_ds_random, delimiter=",")

    print('>>>>>finished plotting data from pdb files!')

    return


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-f", "--folder", help='Folder containing all .pdb files of interest')
    main(vars(parser.parse_args()))

	#print('Rg = {}'.format(Rg(sys.argv[1])))

    # for pdb file in a certain folder
        # calculate Rg values and save in arrays
    # use arrays to plot in a graph with plotly
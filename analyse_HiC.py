import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from transform_gfa import gfa_to_python
import basic_functions as bf   

def short_distance_interactions(fragcontacts, fraglist):
        
    scores = []
    scoreFar = []
    for i in fragcontacts :
        if i[0]+1 == i[1]:
            scores += [i[2]]
        elif  i[1]-i[0] > 10 :
            scoreFar += [i[2]]
    plt.hist(scores, bins=2000)
    plt.xlim ([0,70])
    plt.ylim([0, 1000])
    plt.show()
    
    plt.hist(scoreFar, bins=2000)
    plt.xlim ([0,70])
    plt.ylim([0, 1000])
    plt.show()    
        
def HiC_vs_GFA(hiccontacts, links, fragment_list) :
    
    confirmationOfLinks = [[0 for i in j] for j in links] #a list of list of 0 of the same dimensions as links
    
    for contact in hiccontacts :
        contig1 = fragment_list[contact[0]][0]
        contig2 = fragment_list[contact[1]][0]
        
        for j in range(len(links[contig1*2])):
            if links[contig1*2][j] == contig2*2 or links[contig1*2][j] == contig2*2+1 :
                confirmationOfLinks[contig1*2][j] += contact[2]
                for i in range(len(links[links[contig1*2][j]])):
                    if links[links[contig1*2][j]][i] == contig1*2 :
                        confirmationOfLinks[links[contig1*2][j]][i] += contact[2]                   
                
        for j in range(len(links[contig1*2+1])):
            if links[contig1*2+1][j] == contig2*2 or links[contig1*2+1][j] == contig2*2+1 :
                confirmationOfLinks[contig1*2+1][j] += contact[2]
                for i in range(len(links[links[contig1*2+1][j]])):
                    if links[links[contig1*2+1][j]][i] == contig1*2+1 :
                        confirmationOfLinks[links[contig1*2+1][j]][i] += contact[2]
    
    return confirmationOfLinks

def HiC_vs_GFAtwo(hiccontacts, links, fragment_list, coverage) : #this time we take into account contigs that are two connexions away
    
    confirmationOfLinks = [[0 for i in j] for j in links] #a list of list of 0 of the same dimensions as links
    weightedconfirmationOfLinks = [[0 for i in j] for j in links]
    
    for contact in hiccontacts :
         
        contig1 = fragment_list[contact[0]][0]
        contig2 = fragment_list[contact[1]][0]
        
        for j, neighbor in enumerate(links[contig1*2]):
            #direct neighbor
            if links[contig1*2][j] == contig2*2 or links[contig1*2][j] == contig2*2+1 :
                confirmationOfLinks[contig1*2][j] += contact[2]
                weightedconfirmationOfLinks[contig1*2][j] += contact[2]/coverage[contig1]/coverage[contig2]
                
                for i in range(len(links[links[contig1*2][j]])):
                    if links[links[contig1*2][j]][i] == contig1*2 :
                        confirmationOfLinks[links[contig1*2][j]][i] += contact[2]
                        weightedconfirmationOfLinks[links[contig1*2][j]][i] += contact[2]/coverage[contig1]/coverage[contig2]
            #two connexions away
            for c in range(len(links[neighbor+1-2*neighbor%2])): #we take the other end of the neighbor contig
                if links[neighbor+1-2*neighbor%2][c] == contig2*2 or links[neighbor+1-2*neighbor%2][c] == contig2*2+1 :
                    confirmationOfLinks[contig1*2][j] += contact[2]
                    weightedconfirmationOfLinks[contig1*2][j] += contact[2]/coverage[contig1]/coverage[contig2]
                    
                    confirmationOfLinks[neighbor+1-2*neighbor%2][c] += contact[2]
                    weightedconfirmationOfLinks[neighbor+1-2*neighbor%2][c] += contact[2]/coverage[contig1]/coverage[contig2]
                    
                    for i in range(len(links[links[neighbor+1-2*neighbor%2][c]])) :
                        if links[links[neighbor+1-2*neighbor%2][c]][i] == neighbor+1-2*neighbor%2 :
                            confirmationOfLinks[links[neighbor+1-2*neighbor%2][c]][i] += contact[2]
                            weightedconfirmationOfLinks[links[neighbor+1-2*neighbor%2][c]][i] += contact[2]/coverage[contig1]/coverage[contig2]
                    for i in range(len(links[neighbor])):
                        if links[neighbor][i] == contig1*2 :
                           confirmationOfLinks[neighbor][i] += contact[2]
                           weightedconfirmationOfLinks[neighbor][i] += contact[2]/coverage[contig1]/coverage[contig2]
                    
                
        for j, neighbor in enumerate(links[contig1*2+1]):
            for j in range(len(links[contig1*2+1])):
                if links[contig1*2+1][j] == contig2*2 or links[contig1*2+1][j] == contig2*2+1 :
                    confirmationOfLinks[contig1*2+1][j] += contact[2]
                    weightedconfirmationOfLinks[contig1*2+1][j] += contact[2]/coverage[contig1]/coverage[contig2]
                    
                    for i in range(len(links[links[contig1*2+1][j]])):
                        if links[links[contig1*2+1][j]][i] == contig1*2+1 :
                            confirmationOfLinks[links[contig1*2+1][j]][i] += contact[2]
                            weightedconfirmationOfLinks[links[contig1*2+1][j]][i] += contact[2]/coverage[contig1]/coverage[contig2]
            
            #two connexions away
            otherEndOfNeighbor = neighbor + 1 - 2*(neighbor%2)

            for c in range(len(links[otherEndOfNeighbor])): #we take the other end of the neighbor contig
                if links[otherEndOfNeighbor][c] == contig2*2 or links[otherEndOfNeighbor][c] == contig2*2+1 :
                    confirmationOfLinks[contig1*2+1][j] += contact[2]
                    weightedconfirmationOfLinks[contig1*2+1][j] += contact[2]/coverage[contig1]/coverage[contig2]
                    
                    confirmationOfLinks[otherEndOfNeighbor][c] += contact[2]
                    weightedconfirmationOfLinks[otherEndOfNeighbor][c] += contact[2]/coverage[contig1]/coverage[contig2]

                    for i in range(len(links[links[otherEndOfNeighbor][c]])) :

                        if links[links[otherEndOfNeighbor][c]][i] == otherEndOfNeighbor :
                            confirmationOfLinks[links[otherEndOfNeighbor][c]][i] += contact[2]
                            weightedconfirmationOfLinks[links[otherEndOfNeighbor][c]][i] += contact[2]/coverage[contig1]/coverage[contig2]
                    for i in range(len(links[neighbor])):
                        if links[neighbor][i] == contig1*2 :
                            confirmationOfLinks[neighbor][i] += contact[2]
                            weightedconfirmationOfLinks[neighbor][i] += contact[2]/coverage[contig1]/coverage[contig2]
                    #if verif != 2 :
                     #   print ('il y a un probleme')
    
    return confirmationOfLinks, weightedconfirmationOfLinks

def distance_law(hiccontacts, fragmentList):
    tableDistance = [] # we're going to visualize distance law by looking at the inside of contigs
    tableIntensity = []
    
    for i in hiccontacts :
        if fragmentList[i[0]][0] ==  fragmentList[i[1]][0] and i[2]>0 :#i.e. we are in the same contig
            
            distance = fragmentList[i[1]][1] - fragmentList[i[0]][1]
            tableDistance += [distance]
            tableIntensity += [i[2]]
    
    plt.scatter(np.log10(tableDistance), tableIntensity, alpha = 0.02)
    plt.xlabel('log10(Distance (bp))')
    plt.ylabel('(Intensity)')
    #plt.xlim([0, 150000])
    plt.ylim([0,20])
    plt.show()
    
    #export_to_csv([tableDistance, tableIntensity], 'listsPython/distanceIntensite')
    
def with_how_many_contig_does_one_contig_interact(hiccontactsfile, fragmentList) :
    
    with open(hiccontactsfile) as f:
    
        what_contigs_interact_with_this_one = [[] for i in range(fragmentList[-1][0]+1)]
        for line in f :
            
            line = line.strip('\n')
            line = line.split('\t')
            
            if line[0] != '487796' :#because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]
                
                if contact[2]>4: #this diminishes a lot the number of contig one contig interacts with, probably filtering out errors and rare events
                
                    contig1 = fragmentList[contact[0]][0]
                    contig2 = fragmentList[contact[1]][0]
                
                    if contig2 not in what_contigs_interact_with_this_one[contig1] :
                        what_contigs_interact_with_this_one[contig1] += [contig2]
                    if contig1 not in what_contigs_interact_with_this_one[contig2] :
                        what_contigs_interact_with_this_one[contig2] += [contig1]
        
        how_many_contigs_interact_with_this_one = [len(x) for x in what_contigs_interact_with_this_one]
        plt.hist(how_many_contigs_interact_with_this_one, bins = 20)
        #plt.xlim([0,200])
        plt.xlabel('Number of other contigs interacting with one at least by 4 contacts')
        plt.ylabel('Number of contig interacting with x others')
        plt.show()

#here comes the neutral test for our test_HiC_vs_GFA : we're going to break down contigs and see how much HiC contact they have
def testHiC_vs_GFA(hiccontacts, info_contigs) :
    
    contactNumber = 0
    score = []
    
    for contig in range(len(info_contigs)):
        
        start_frag = info_contigs[contig][3]
        end_frag = info_contigs[contig][3] + info_contigs[contig][2]
        
        if info_contigs[contig][2] > 29 :
            cut = random.randint(15, info_contigs[contig][2]-14) #we cut the contig in two random parts
        
            score += [0]
            
            if contig < 5 :
                print(contig, start_frag, end_frag, cut)
            
            while hiccontacts[contactNumber][0] < start_frag + cut :
                
                if hiccontacts[contactNumber][0] > start_frag and\
                    hiccontacts[contactNumber][1] >= start_frag+cut and hiccontacts[contactNumber][1] < end_frag :
                        score[-1] += hiccontacts[contactNumber][2]
                
                contactNumber += 1
        
        while hiccontacts[contactNumber][0] < end_frag and contactNumber < len(hiccontacts)-1:
            contactNumber += 1
    
    print(score)
    plt.hist(score)
    
def detect_fishy_links(links, confirmationOfLinks, coverage):
    
    #we're going to detect, when there is an ambiguity, if one path looks unlikely
    badlinks = []*len(links)
    for i in coverage : #to ensure we don't get absurdly high multiplicative factors
        if i < 0.01 :
            i = 0.01

    for endOfContig in range(len(links)) :
                
        if len(links[endOfContig]) > 1 :
            
            weightedConfirmation = [-1]*len(links[endOfContig])
            for i in range(len(links[endOfContig])) :
                
                if coverage[int(endOfContig/2)]>0.01 and coverage[int(links[endOfContig][i]/2)]>0.01 :
                    weightedConfirmation[i] = confirmationOfLinks[endOfContig][i]/coverage[int(links[endOfContig][i]/2)]
                    
            maximum = np.max(weightedConfirmation)
            reliable = -1 not in weightedConfirmation
        
            if reliable :
                for i in range(len(weightedConfirmation)) :
                    if weightedConfirmation[i]<0.1*maximum :
                        
                        badlinks += [[endOfContig, i]]
                        print('There is a suspect link here : ', endOfContig, links[endOfContig], \
                              confirmationOfLinks[endOfContig],\
                              weightedConfirmation,)
    return badlinks
     
def interactionMatrix(hiccontactsfile, fragmentList, coverage):
    with open(hiccontactsfile) as f:
    
        interactionMatrix = [[0 for i in range(fragmentList[-1][0]+1)] for j in range(fragmentList[-1][0]+1)]
        for line in f :
            
            line = line.strip('\n')
            line = line.split('\t')
            
            if line[0] != '487796' :#because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]
                                
                contig1 = fragmentList[contact[0]][0]
                contig2 = fragmentList[contact[1]][0]
                
                interactionMatrix[contig1][contig2] += contact[2]*1000000/coverage[contig1]/coverage[contig2] #*1000000 so that the numbers are not too small (the risk being their being considered 0)
                interactionMatrix[contig2][contig1] += contact[2]*1000000/coverage[contig1]/coverage[contig2]
    return interactionMatrix

#hiccontacts = read_abs_fragments_contact_weighted('data/results/abs_fragments_contacts_weighted.txt')
#hiccontacts = import_from_csv('listsPython/hiccontacts.csv')
#print(hiccontacts[:20])
#print(hiccontacts[:100])
fragmentList = bf.read_fragment_list('data/results/fragments_list.txt')
#print(fragmentList[:100])
#infcontigs = read_info_contig('data/results/info_contigs.txt')
#links = gfa_to_python(1312)
#short_distance_interactions(hiccontacts, fragmentList)

#links = import_from_csv('listsPython/links.csv')
#links = [[int(i) for i in j] for j in links]    

#hiccontacts = import_from_csv('listsPython/hiccontacts.csv')

#confirmationOfLinks = HiC_vs_GFAtwo('data/results/abs_fragments_contacts_weighted.txt', links, fragmentList)

#distance_law(hiccontacts, fragmentList)
#testHiC_vs_GFA(hiccontacts, infcontigs)
#determine_HiC_coverage(hiccontacts, infcontigs, fragmentList)
#confirmationOfLinks = import_from_csv('listsPython/confirmationsDeslinks.csv')
#print(confirmationOfLinks[:19])

#check_links(links)
#coverage = determine_HiC_coverage(hiccontacts, infcontigs, fragmentList)
coverage = bf.import_from_csv('listsPython/HiCcoverage.csv')
coverage = [x[0] for x in coverage]

#conf, confweight = HiC_vs_GFAtwo('data/results/abs_fragments_contacts_weighted.txt', links, fragmentList, coverage)
#print(conf[:20],confweight[:20])
#print(links[:20])

#with_how_many_contig_does_one_contig_interact('data/results/abs_fragments_contacts_weighted.txt', fragmentList)

interaction_Matrix = interactionMatrix('data/results/abs_fragments_contacts_weighted.txt', fragmentList, coverage)
bf.export_to_csv(interaction_Matrix, 'listsPython/interactionMatrix.csv')
print(interaction_Matrix[217][323], interaction_Matrix[217][359])

print('Finished')

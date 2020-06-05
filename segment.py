#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.integrate as integrate

class Segment:
    def __init__(self, segIndex, segListOfContig, segNamesOfContig, segOrientationOfContigs, segLengths, segInsideCIGARs = None, segLinks = [[],[]], segOtherEndOfLinks = [[],[]], lock = False):
        
        if len(segLinks[0]) != len(segOtherEndOfLinks[0]) or len(segLinks[1]) != len(segOtherEndOfLinks[1]) :
            print('ERROR in the links while initializing a segment')
            return 0
        
        if any(i!=0 and i!=1 for i in segOtherEndOfLinks[0]) or any(i!=0 and i!=1 for i in segOtherEndOfLinks[1]):
            print('ERROR in the links while initializing a segment')
            return 0
        
        if len(segListOfContig) != len(segOrientationOfContigs) :
            print('ERROR in initializing the orientations of contigs within a segment')
            return 0
        
        if segInsideCIGARs == None :
            segInsideCIGARs = ['*' for i in range(len(segListOfContig)-1)]
            
        self._index = segIndex
        
        #this group of attributes are linked array : element n in one corresponds with element n in the other. Therefore they shouldn't be modified independantly
        self._namesOfContigs = segNamesOfContig.copy() #names are strings with which sequences are described in the GFA
        self._listOfContigs = segListOfContig.copy() #listOfContigs are numbers with which contigs are referenced in interactioMatrix
        self._orientationOfContigs = segOrientationOfContigs.copy() #1 being '+' orientation, 0 the '-' orientation
        self._lengths = segLengths.copy()
        self._insideCIGARs = segInsideCIGARs.copy()
        
        #this group of attribute are linked arrays : one should never be modified without the others 
        self._links = [[],[]] #two lists of segments with which the segment is linked, at the left end and at the right end
        self._otherEndOfLinks = [[],[]] #for each link, indicates the side of the other segment on which the link arrives
        self._CIGARs = [[],[]]
        
        self._freezed = [False, False]
        self._locked = lock #That is to duplicate a contig only once in each meerge_contigs
        
    # getters
    
    def get_index(self):
        return self._index

    def get_listOfContigs(self):
        return self._listOfContigs
    
    def get_orientations(self):
        return self._orientationOfContigs
    
    def get_insideCIGARs(self):
        return self._insideCIGARs
    
    def get_links(self):
        return self._links
    
    def get_otherEndOfLinks(self):
        return self._otherEndOfLinks
    
    def get_CIGARs(self):
        return self._CIGARs
    
    def get_lengths(self):
        return self._lengths
    
    def get_namesOfContigs(self):
        return self._namesOfContigs
    
    def get_freezed(self):
        return self._freezed
    
    def get_locked(self):
        return self._locked
    
    def print_complete(self):
        print(self._namesOfContigs, [s.names for s in self._links[0]], [s.names for s in self._links[1]])
    
    # setters 
    
    def set_index(self, newIndex):
        self._index = newIndex
    
    def set_listOfContigs(self, newListOfContigs):
        self._listOfContigs = newListOfContigs  
        
    def freeze(self, endOfSegment): 
        self._freezed[endOfSegment] = True
  
    def freezeNode(self, endOfSegment):
        self._freezed[endOfSegment] = True
        for n, neighbor in enumerate(self._links[endOfSegment]) :
            neighbor.freeze(self._otherEndOfLinks[n])
        
    def unfreeze(self):
        self._freezed = [False, False]
        
    def set_locked(self, b):
        self._locked = b
        
    def hash(self) : #a fucntion assigning an int to a segment, for example useful in export_to_GFA
            h = int(''.join([str(i) for i in self._listOfContigs]))
            return h
    # properties
    
    index = property(get_index, set_index)
    
    names = property(get_namesOfContigs)
    listOfContigs = property(get_listOfContigs, set_listOfContigs)
    orientations = property(get_orientations)
    lengths = property(get_lengths)
    insideCIGARs = property(get_insideCIGARs)
    
    links = property(get_links)
    otherEndOfLinks = property(get_otherEndOfLinks)
    CIGARs = property(get_CIGARs)
    
    freezed = property(get_freezed)
    locked = property(get_locked, set_locked)

    #other functions that handle segments
    
    #function which goals is to return the intensity of HiC contacts between another segment and this one
    def interaction_with_contigs(self, segment, interactionMatrix, dist_law = lambda x:1, copiesnumber = None, commonContigs = [], bestSignature = 1000):
        
        if copiesnumber == None :
            copiesnumber = [1 for i in interactionMatrix]
            
        absoluteScore = 0
        relativeScore = 0
        
        partial_area = 0
        orientation = -1 # if supercontig is directly linked to the candidates, then this variable tells us by which end
        
        if segment in self.links[0]:
            orientation = 0
        elif segment in self.links[1]:
            orientation = 1
        else: #if the supercontigs are not touching, computing the partial area is useless, but harmless
            print('ERROR : trying to compute an interaction with supercontigsaretouching=True but actually not True')
            
        # lengthForNow is a variable keeping track of how far the examined contig of the listOfSuperContigs is from supercontig
        lengthForNow = orientation * np.sum([i for i in self.lengths])
        
        for contig in range(len(self.listOfContigs)) :

            newLengthForNow = (lengthForNow + self.lengths[contig] * (0.5 - orientation) * 2)
            if self.listOfContigs[contig] not in commonContigs:
                # computing the partial area : that way, small supercontigs are not penalized when compared to much longer ones
                partial_area += np.abs(integrate.quad(dist_law, newLengthForNow, lengthForNow)[0])

            for contigInSegment in segment.listOfContigs:
                if self.listOfContigs[contig] not in commonContigs and copiesnumber[contigInSegment] <= bestSignature:
                    absoluteScore += interactionMatrix[contigInSegment][self.listOfContigs[contig]]
                    relativeScore += interactionMatrix[contigInSegment][self.listOfContigs[contig]]
                else:
                    absoluteScore += interactionMatrix[contigInSegment][self.listOfContigs[contig]]

            lengthForNow = newLengthForNow
            
            return absoluteScore, relativeScore, partial_area
    
    def add_link_from_GFA(self, GFAline, names, segments, leftOrRight) : #leftOrRight = 0 when the segment is at the beginning of a link (left of a GFA line), 1 otherwise
        
        l = GFAline.strip('\n').split('\t')
 
        if len(l) < 5 :
            print('ERROR : expected at least 5 fields in line ', GFAline)
        
        if l[0] != 'L':
            print('ERROR : trying to add a link from a GFA line that does not start with "L"')
        
        else :
            o1,o2 = -1, -1
            
            if l[2] == '-':
                o1 = 0
            elif l[2] == '+' :
                o1 = 1
                
            if l[4] == '-':
                o2 = 0
            elif l[4] == '+' :
                o2 = 1
                
            if o1 == -1 or o2 == -1 :
                print('ERROR while creating a link : orientations not properly given.')
                print('Problematic line : ', GFAline)
                
            if leftOrRight == 0 and o1 != self._orientationOfContigs[0] :
                self._links[0].append(segments[names.index(l[3])])
                self._otherEndOfLinks[0].append(1-o2)
                if len(l) > 5 :
                    self._CIGARs[0].append(l[5])
                else :
                    self._CIGARs[0].append('*')
                    
            elif leftOrRight == 0 and o1 == self._orientationOfContigs[-1] :
                self._links[1].append(segments[names.index(l[3])])
                self._otherEndOfLinks[1].append(1-o2)
                if len(l) > 5 :
                    self._CIGARs[1].append(l[5])
                else :
                    self._CIGARs[1].append('*')
                
            elif leftOrRight == 1 and o2 == self._orientationOfContigs[0] :
                self._links[0].append(segments[names.index(l[1])])
                self._otherEndOfLinks[0].append(o1)
                if len(l) > 5 :
                    self._CIGARs[0].append(l[5])
                else :
                    self._CIGARs[0].append('*')
                    
            elif leftOrRight == 1 and o2 != self._orientationOfContigs[-1] :
                self._links[1].append(segments[names.index(l[1])])
                self._otherEndOfLinks[1].append(o1)
                if len(l) > 5 :
                    self._CIGARs[1].append(l[5])
                else :
                    self._CIGARs[1].append('*')
                    
            else :
                print('ERROR while trying to add a new link from the gfa : could not locate a correct name')
                         
    #this adds the end of a links, but only on this segment, not on the other end
    def add_end_of_link(self, endOfSegment, segment2, endOfSegment2, CIGAR = '*'):
        self._links[endOfSegment].append(segment2)
        self._otherEndOfLinks[endOfSegment].append(endOfSegment2)
        self._CIGARs[endOfSegment].append(CIGAR.copy())

    def remove_end_of_link(self, endOfSegment, segmentToRemove, endOfSegmentToRemove = None): #endOfSegmentToRemove is there in case there exists two links between self[endOfSegment] and segment to remove. Needed for extra security
        
        #first determine the index of the segment to remove
        index = self._links[endOfSegment].index(segmentToRemove)

        if endOfSegmentToRemove != None :
            if self._links[endOfSegment].count(segmentToRemove) == 2 :   
                for i in range(len(self._links[endOfSegment])):
                    if self._links[endOfSegment][i] == segmentToRemove and self._otherEndOfLinks[endOfSegment][i] == endOfSegmentToRemove :
                        index = i
            
        #then remove the end of unwanted link in all attributes
        del self._links[endOfSegment][index]
        del self._otherEndOfLinks[endOfSegment][index]
        del self._CIGARs[endOfSegment][index]
            
#This function is OUTSIDE the class. It takes two segments and the end of the first segment which is linked to the second. It appends a merged contig to the listOfSegments, without modifying the two inputed segments
def merge_two_segments(segment1, endOfSegment1, segment2, listOfSegments):
    
    if segment1.links[endOfSegment1].count(segment2) > 1 : #this means a loop
        return 0
      
    # creating a new segment
    orientation1 = endOfSegment1*2-1
    endOfSegment2 = segment1.otherEndOfLinks[endOfSegment1][segment1.links[endOfSegment1].index(segment2)]
    orientation2 = -2*endOfSegment2+1
    
    orientationOfContigs1 = segment1.orientations
    orientationOfContigs2 = segment2.orientations
    
    if orientation1 == -1 : #then change the orientation of all the contigs within the segment, since the segment will be mirrored in the new supersegment
        orientationOfContigs1 = [1-i for i in orientationOfContigs1]
    if orientation2 == -1 :
        orientationOfContigs2 = [1-i for i in orientationOfContigs2]
        
    CIGAR = segment1.CIGARs[endOfSegment1][segment1.links[endOfSegment1].index(segment2)]
    
    newSegment = Segment(len(listOfSegments), \
                         segment1.listOfContigs[::orientation1] + segment2.listOfContigs[::orientation2],\
                             segment1.names[::orientation1] + segment2.names[::orientation2],\
                                orientationOfContigs1+orientationOfContigs2,\
                                segment1.lengths[::orientation1]+segment2.lengths[::orientation2], \
                                segment1.insideCIGARs + [CIGAR] + segment2.insideCIGARs,\
                                [segment1.links[1-endOfSegment1], \
                                segment2.links[1-endOfSegment2]], \
                                [segment1.otherEndOfLinks[1-endOfSegment1], \
                                segment2.otherEndOfLinks[1-endOfSegment2]],\
                                lock = True)
    
    listOfSegments.append(newSegment)
    
    #building the other end of links with the new segment
    for n, neighbor in enumerate(newSegment.links[0]) :
        neighbor.add_end_of_link(neighbor.otherEndOfLinks[0][n], newSegment, 0, CIGAR = newSegment.CIGARs[0][n])
    for n, neighbor in enumerate(newSegment.links[1]) :
        neighbor.add_end_of_link(neighbor.otherEndOfLinks[1][n], newSegment, 1, CIGAR = newSegment.CIGARs[1][n])

#function creating a link between two ends of contigs, OUTSIDE of the class
def add_link(segment1, end1, segment2, end2, CIGAR = '*'):
    segment1.add_end_of_link(end1, segment2, end2, CIGAR)
    segment2.add_end_of_link(end2, segment1, end1, CIGAR)
           


## A few lines to test the functions of the file

# s1 = Segment(0, [0,1], [1,1], [1000])
# s2 = Segment(1, [2], [1], [500])
# add_link(s1,1,s2,0,'60M')

# listOfSegments = [s1, s2]

# merge_two_segments(s1, 1, s2, listOfSegments)

# print([i.listOfContigs for i in listOfSegments])
# print([i.orientations for i in listOfSegments])

                    
                    
                    
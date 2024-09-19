#recursively see if we can reach 100000 from here
def length_from_here(segment, end, length_now, max_length) :
        
        if length_now > max_length :
            return length_now
        
        if len(segment.links[end]) == 0 :
            return length_now + segment.length
        
        best_length = 0
        for n, neighbor in enumerate(segment.links[end]) :
            length_here = length_from_here(neighbor, 1-segment.otherEndOfLinks[end][n], length_now + segment.length, max_length)
            if length_here > best_length :
                best_length = length_here
            if length_here > max_length :
                return length_here
            
        return best_length

#in case the graph is noisy, determine haploid coverage with long contigs only, then pop the bubbles and delete the tips     
def clean_graph(segments, genome_size) :

    coverage_total = 0
    for s in segments :
        coverage_total += s.length * s.depths[0]

    haploid_coverage = coverage_total / genome_size
        
    #now pop the bubbles and trim the tips
    contigs_to_delete = []
    for s in segments:

        #check if the contig is in a dead end
        if len(s.links[0]) == 0 and len(s.links[1]) == 1:
            #is this contig not the best option ? GO through the neighbors of the neighbor
            neighbor = s.links[1][0]
            best_coverage = 0
            best_length = length_from_here(s.links[1][0], s.otherEndOfLinks[1][0], s.length, 100000)
            for n, neighbor2 in enumerate(s.links[1][0].links[s.otherEndOfLinks[1][0]]) :
                if neighbor2.ID != s.ID :
                    if neighbor2.depths[0] > best_coverage :
                        best_coverage = neighbor2.depths[0]
                        best_length = length_from_here(neighbor2, 1-s.otherEndOfLinks[1][0], s.length, 100000)
            
            if s.length < 0.2*best_length and s.depths[0] < best_coverage and s.depths[0] < haploid_coverage*0.5 :
                contigs_to_delete.append(s.ID)
                
        if len(s.links[1]) == 0 and len(s.links[0]) == 1:
            #is this contig not the best option ? GO through the neighbors of the neighbor
            neighbor = s.links[0][0]
            best_coverage = 0
            best_length = length_from_here(s.links[0][0], s.otherEndOfLinks[0][0], s.length, 100000)
            for n, neighbor2 in enumerate(s.links[0][0].links[s.otherEndOfLinks[0][0]]) :
                if neighbor2.ID != s.ID :
                    if neighbor2.depths[0] > best_coverage :
                        best_coverage = neighbor2.depths[0]
                        best_length = length_from_here(neighbor2, 1-s.otherEndOfLinks[0][0], s.length, 100000)
            
            if s.length < 0.2*best_length and s.depths[0] < best_coverage and s.depths[0] < haploid_coverage*0.5 :
                contigs_to_delete.append(s.ID)

        #check if it is a bubble: one neighbor left and right, this neighbor has only one other neighbor, which is the same for both
        if len(s.links[0]) == 1 and len(s.links[1]) == 1:
            left_neighbor = s.links[0][0]
            left_end = s.otherEndOfLinks[0][0]
            right_neighbor = s.links[1][0]
            right_end = s.otherEndOfLinks[1][0]
            if len(left_neighbor.links[left_end]) == 2 and len(right_neighbor.links[right_end]) == 2:
                other_neighbor_left = left_neighbor.links[left_end][0]
                if other_neighbor_left.ID == s.ID :
                    other_neighbor_left = left_neighbor.links[left_end][1]
                other_neighbor_right = right_neighbor.links[right_end][0]
                if other_neighbor_right.ID == s.ID :
                    other_neighbor_right = right_neighbor.links[right_end][1]

                if other_neighbor_left.ID == other_neighbor_right.ID : #then this is a bubble, let us see if we should pop it
                    #check if the bubble is an area of 
                    if other_neighbor_left.depths[0] > s.depths[0] and s.depths[0] < haploid_coverage*0.6  \
                        and right_neighbor.depths[0] < haploid_coverage*1.5 \
                        and left_neighbor.depths[0] < haploid_coverage*1.5 :
                        contigs_to_delete.append(s.ID)
                    elif other_neighbor_left.depths[0] > s.depths[0] and s.depths[0] < haploid_coverage*0.2 :
                        contigs_to_delete.append(s.ID)


    #now delete the contigs
    toDetete_set = set(contigs_to_delete)
    toDelete_idx = []
    idx = 0
    for s in segments :
        if s.ID in toDetete_set :
            s.cut_all_links()
            toDelete_idx.append(idx)
        idx += 1

    for i in toDelete_idx[::-1] :
        del segments[i]
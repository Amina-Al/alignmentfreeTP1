from loading import load_directory
from kmers import stream_kmers, kmer2str
import numpy as np



def jaccard(fileA, fileB, k):
     # dico kmers
    dico_kmer = {}
    intersection = 0
    union = 0
    
    # Get the ones in seqA
    for one_kmer in stream_kmers(fileA,k):
        # Put in dico if not already in it
        if one_kmer not in dico_kmer:
            union += 1
            dico_kmer[one_kmer] = 1
        else:
            union += 1
            dico_kmer[one_kmer] +=1
    for one_kmer in stream_kmers(fileB,k):
        # Put in dico if not already in it
        if one_kmer in dico_kmer:
            intersection +=1
            dico_kmer[one_kmer] -= 1
            if dico_kmer[one_kmer] == 0:
                del dico_kmer[one_kmer]
        else:
            union += 1   
             
    return intersection/union


if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21 
    dist_matrix = np.identity((len(files)))

    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            coef = jaccard(files[filenames[i]],files[filenames[j]],k)
            dist_matrix[i,j] = coef
            #print(filenames[i], filenames[j], coef)

    print(dist_matrix)
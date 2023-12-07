
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    """
    Input : text : from which we want the kmers
            k : length of the kmer 
    Output : kmer"
    """
    
    # Dictionnary to encode and have the complementary kmers
    encodage = {'A':0b00,'C':0b01,'T':0b10,'G':0b11}
    complementary = {'A':'T','C':'G', 'T':'A', 'G':'C'}
    
    #Initialisation kmerb 
    kmerb = 0b00
    compl_kmerb = 0b00
    
    # Initialisation of a list that contain all the kmer found 
    list_kmers = []
    #print(f"IL ME SOULE CE TRUC 1 {type(kmerb)}")
    # We get the k first values of the text 
    for let in text[0:k]:
        #print(f"IL ME SOULE CE TRUC 2 {type(kmerb)}")
        kmerb = kmerb << 2
        compl_kmerb <<= 2
        #print(f"IL ME SOULE CE TRUC 3 {type(kmerb)}")
        kmerb= kmerb | encodage[let]
        compl_kmerb |= encodage[complementary[let]]
        #print(f"IL ME SOULE CE TRUC 4 {type(kmerb)}")
    
    print(f"kmerb: {kmerb} et son complementary: {compl_kmerb}")
    # Choose the minimum between 2 complementary kmer
    min_kmerb = min(kmerb,compl_kmerb)
    list_kmers.append(min_kmerb)
    yield min_kmerb
    
    #Calcul mask
    mask = (1<<(k*2))-1
    
    
    for nucl in range(k,len(text)):
        kmerb <<= 2
        compl_kmerb <<= 2
        kmerb |= encodage[text[nucl]]
        compl_kmerb |= encodage[complementary[text[nucl]]]
        kmerb &= mask
        compl_kmerb &= mask
        if min(kmerb,compl_kmerb) not in list_kmers:
            yield min(kmerb,compl_kmerb)
    
    
    


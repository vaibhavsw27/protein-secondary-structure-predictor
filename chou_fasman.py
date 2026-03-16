# Chou Fasman Secondary Structure Prediction Code
# By : Vaibhav Samrat Waghaye 2024597



# Residue propensities (Pα, Pβ)

chou_fasman ={
    "A": [1.45, 0.97], #Ala
    "R": [0.79, 0.90], #Arg
    "N": [0.73, 0.65], #Asn
    "D": [0.98, 0.80], #Asp
    "C": [0.77, 1.30], #Cys
    "Q": [1.17, 1.23], #Gln
    "E": [1.53, 0.26], #Glu
    "G": [0.53, 0.81], #Gly
    "H": [1.24, 0.71], #His
    "I": [1.00, 1.60], #Ile
    "L": [1.34, 1.22], #Leu
    "K": [1.07, 1.19], #Lys
    "M": [1.20, 1.67], #Met
    "F": [1.12, 1.28], #Phe
    "P": [0.59, 0.62], #Pro
    "S": [0.79, 0.72], #Ser
    "T": [0.82, 1.20], #Thr
    "W": [1.14, 1.19], #Trp
    "Y": [0.61, 1.29], #Tyr
    "V": [1.14, 1.65] #Val
}






def helical_regions(seq,size):
    helix_regions=["-" for _ in range(size)] #intializing helix regions list


    for i in range(size-5):
        window = seq[i:i+6] #6 residue window

        count = 0 #count of residues with propensity >=1.00
        for x in window:
            if chou_fasman.get(x)[0]>=1.00: #checking propensity has to be >=1.00 
                #Note- values like >=1.03 can also be used based on sources like wikipedia,national library of medicine etc.
                # I have used >=1.00 as per class notes , but either can be used
    
                count += 1

        if count >= 4:

            start=i
            end=i+6

            #extending to the right
            while end<size:
                window = seq[end-3:end+1]
                total = 0
                for x in window:
                    total += chou_fasman[x][0]
                avg = total / 4
                if avg >= 1.00:
                    end += 1
                else:
                    break

            #extending to the left
            while start > 0:
                window = seq[start-1:start+3]
                total = 0
                for x in window:
                    total += chou_fasman[x][0]
                avg = total / 4
                if avg >= 1.00:
                    start -= 1
                else:
                    break

            for j in range(start,end):
                helix_regions[j] = "H" #marking helix regions


    return helix_regions



def beta_regions(seq, size):
    strand_regions = ["-" for _ in range(size)] #initializing beta strand regions list

    for i in range(size - 4):
        window = seq[i:i+5] #5 residue window

        count = 0  #count of residues with propensity >=1.00
        for x in window:
            if chou_fasman.get(x)[1] >= 1.00: #checking propensity has to be >=1.00
                count += 1

        if count >= 3:  # nucleation site
            start = i
            end = i + 5

            # extend to the right
            while end < size:
                window = seq[end-3:end+1]
                total = 0
                for x in window:
                    total += chou_fasman[x][1]
                avg = total / 4
                if avg >= 1.00:
                    end += 1
                else:
                    break

            # extend to the left
            while start > 0:
                window = seq[start-1:start+3]
                total = 0
                for x in window:
                    total += chou_fasman[x][1]
                avg = total / 4
                if avg >= 1.00:
                    start -= 1
                else:
                    break

            for j in range(start, end):
                strand_regions[j] = "S" #marking beta strand regions

    return strand_regions




def resolve_conflicts(seq, H, S,size):
    
    final = [] #final secondary structure list

    for i in range (size):
        x=seq[i]

        if H[i] == "H" and S[i] == "S":
            if chou_fasman[x][0]>chou_fasman[x][1]: #resolving conflicts based on higher propensity
                final.append("H")
            else:
                final.append("S")

        elif H[i] == "H":
            final.append("H")
        elif S[i] == "S":
            final.append("S")
        else:
            final.append("-") 

    return final



def get_ranges(structure,letter,seq,size):
    regions = [] #list to store regions
    i = 0
    while i<size:
        if structure[i]==letter: #checking whether we found the start of a region
            start=i
            while i<size and structure[i]==letter:
                i += 1

            end=i #end of the region
            residues = seq[start:end]
            regions.append((start + 1,end, residues)) #storing 1-based index
        else:
            i += 1
    return regions



def main():

    seq=input("Please enter your sequence- ")
    #below is an example sequence , which you can use
    # seq = """MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYLEKPMEIARIVARCLWEESRLLQTAATAAQQGGQANHPTAAVVTEKQQMLEQHLQDVRKRVQDLEQKMKVVENLQDDFDFNYKTLKSQGDMQDLNGNNQSVTRQKMQQLEQMLTALDQMRRSIVSELAGLLSAMEYVQKTLTDEELADWKRRQQIACIGGPPNICLDRLENWITSLAESQLQTRQQIKKLEELQQKVSYKGDPIVQHRPMLEERIVELFRNLMKSAFVVERQPCMPMHPDRPLVIKTGVQFTTKVRLLVKFPELNYQLKIKVCIDKDSGDVAALRGSRKFNILGTNTKVMNMEESNNGSLSAEFKHLTLREQRCGNGGRANCDASLIVTEELHLITFETEVYHQGLKIDLETHSLPVVVISNICQMPNAWASILWYNMLTNNPKNVNFFTKPPIGTWDQVAEVLSWQFSSTTKRGLSIEQLTTLAEKLLGPGVNYSGCQITWAKFCKENMAGKGFSFWVWLDNIIDLVKKYILALWNEGYIMGFISKERERAILSTKPPGTFLLRFSESSKEGGVTFTWVEKDISGKTQIQSVEPYTKQQLNNMSFAEIIMGYKIMDATNILVSPLVYLYPDIPKEEAFGKYCRPESQEHPEADPGSAAPYLKTKFICVTPTTCSNTIDLPMSPRTLDSLMQFGNNGEGAEPSAGGQFESLTFDMELTSECATSPM"""
    
    size=len(seq) #size of sequence

    # raw predictions
    raw_H = helical_regions(seq,size)
    raw_S = beta_regions(seq,size)



    #rfinal predictions after resolving conflicts
    final = resolve_conflicts(seq,raw_H,raw_S,size)


    # final derived maps
    final_H = []
    for conflict in final:
        if conflict == "H":
            final_H.append("H")
        else:
            final_H.append("-")

    final_S = []
    for c in final:
        if c == "S":
            final_S.append("S")
        else:
            final_S.append("-")



    #getting regions
    raw_H_regions=get_ranges(raw_H, "H", seq,size)
    final_H_regions=get_ranges(final_H, "H", seq,size)
    raw_S_regions=get_ranges(raw_S, "S", seq,size)
    final_S_regions=get_ranges(final_S, "S", seq,size)




    #printing the outputs
    print("a)\n")
    print("i) Raw Helical Regions:")
    for s, e, residues in raw_H_regions:
        print(f"{s}-{e}: {residues}")


    print("\nii) Final Helical Regions:")
    for s, e, residues in final_H_regions:
        print(f"{s}-{e}: {residues}")


    print("\niii) Raw Helical Map:")
    print("".join(raw_H))
    print("\niv): Final Helical Map:")
    print("".join(final_H))

    
    print("\nb)")

    print("i) Raw Beta Strand Regions:")
    for s, e, residues in raw_S_regions:
        print(f"{s}-{e}: {residues}")    

    print("\nii) Final Beta Strand Regions:")
    for s, e, residues in final_S_regions:
        print(f"{s}-{e}: {residues}")
    
    print("\niii) Raw Beta Strand Map:")
    print("".join(raw_S))

    print("\niv) Final Beta Strand Map:")
    print("".join(final_S))


    print("\nc)")
    print("\nFinal Secondary Structure:")
    print("".join(final))


main()
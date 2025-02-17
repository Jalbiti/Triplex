from Bio import SeqIO
import math

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def calculate_tm_GP(seq, pH = 7.0, ct = 12.0, cd = 0.05, type_TFO = 'RNA', mismatch = 0, penalty = True):
  R = 0.0019858
  occ_a = occurrences(seq,"A")
  # Strictly accept triplex triads
  if (occ_a > mismatch):
      return(0)
  else:  
      occ_c = occurrences(seq,"C")
      occ_cc = occurrences(seq,"CC")

      if (type_TFO == 'RNA'):
          seq = seq.replace('T', 'U')

          occ_u = occurrences(seq,"U")
          occ_uu = occurrences(seq,"UU")
          occ_uc = occurrences(seq,"UC")
          occ_cu = occurrences(seq,"CU")

          # RNA
          enthalpy = -10.95*(occ_cc)-5.73*(occ_uc+occ_cu)-6.44*(occ_uu)
          free_energy = -1.891*(occ_c)-0.758*(occ_u)-0.331*(occ_cc)+2.646+(occ_c)*(pH-5.6)*(0.893-0.005*(occ_cc))

      elif (type_TFO == 'DNA'):
          occ_t = occurrences(seq,"T")
          occ_tt = occurrences(seq,"TT")
          occ_tc = occurrences(seq,"TC")
          occ_ct = occurrences(seq,"CT")

          # DNA
          enthalpy = -4.9*(occ_cc)-8.9*(occ_tc+occ_ct)-7.4*(occ_tt)
          free_energy = -3.00*(occ_c)-0.65*(occ_t)+1.65*(occ_cc)+6.0+(occ_c)*(pH-5.0)*(1.26-0.08*(occ_cc))
  
      # Tm
      # Compute the Tm based on DG, DH and the concentration
      # tm = (298 * DH / (DH-DG - 298*R*log(4/Ct)))
      # Ct: duplex + third strand concentration
      ref_t = 310
      tm_trip = (ref_t * enthalpy / (enthalpy - free_energy - ref_t * R * math.log(4/ct/1e-6)) )
      tm_trip_c = (tm_trip - 273.15)

      # c50
      # Compute concentration of TFO needed to achieve 50% of triplex formation.
      # o_c50 = (1/K_eq) + d/2 = exp(DG/kT) + d/2
      ref_t = 310
      oligo_c50 = 1e6 * math.exp(free_energy / (ref_t * R)) + 0.5 * cd

      if(penalty):
        occ_g = occurrences(seq,"G")
        ngaps = occ_a + occ_g
        tm_trip_c = tm_trip_c - 10*ngaps

      return(tm_trip_c)

# Define the minimum and maximum length for DNA sequences
m = 10  # Minimum length
l = 30  # Maximum length
desired_Tm = 30
#desired_Tm = 45 # For very stable TFOs

# Open the FASTA file and extract sequences with desired Tm

fasta_file = "data/miRNA.fa"
outfile = 'miRNA.out'

with open(fasta_file, "r") as handle:
    records = SeqIO.parse(fasta_file, "fasta")
    fout = open(outfile, "w")
    for record in records: 
        sequence = str(record.seq)
        # save id to map back
        id = str(record.id)
        seq_length = len(sequence)
        # Extract all possible sequences from the read
        for start in range(seq_length):
            for end in range(start + m, min(seq_length, start + l) + 1):
                subsequence = sequence[start:end]
                Tm = calculate_tm_GP(subsequence)
                # Check if the subsequence meets the Tm criteria
                if Tm > desired_Tm:
                    #fout.write("{}\n".format(subsequence))
                    fout.write("{}\t{}\n".format(subsequence, id))
    fout.close()
handle.close()


fasta_file = "data/gencode.v41.lncRNA_transcripts.fa"
outfile = 'lncRNA.out'

with open(fasta_file, "r") as handle:
    records = SeqIO.parse(fasta_file, "fasta")
    fout = open(outfile, "w")
    for record in records:
        sequence = str(record.seq)
        seq_length = len(sequence)
        # save id to map back
        id = str(record.id)
        # Extract all possible sequences from the read
        for start in range(seq_length):
            for end in range(start + m, min(seq_length, start + l) + 1):
                subsequence = sequence[start:end]
                Tm = calculate_tm_GP(subsequence)
                # Check if the subsequence meets the Tm criteria
                if Tm > desired_Tm:
                    #fout.write("{}\n".format(subsequence))
                    fout.write("{}\t{}\n".format(subsequence, id))

    fout.close()
handle.close()


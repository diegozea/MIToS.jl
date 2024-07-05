using MIToS.MSA
using MIToS.Utils

# Define a new file format type for A3M files
struct A3M <: FileFormat end

# Function to handle inserts in sequences, modifying sequences to include insertions
function _add_inserts(SEQS)
    seq_len = length.(SEQS)  
    ncol = maximum(seq_len)  
    j = 1
    is_insert = false  
    
    # Iterate over columns of the sequence alignment
    while j <= ncol
        for i in 1:length(SEQS)
            if j <= seq_len[i]
                res = SEQS[i][j]  
                if islowercase(res) || res == '.'
                    is_insert = true  
                    break
                end
            end
        end
        
        if is_insert
            for i in 1:length(SEQS)
                seq = SEQS[i]
                if j > seq_len[i]
                    SEQS[i] = seq * "." 
                    seq_len[i] += 1
                else
                    res = seq[j]
                    if isuppercase(res) || res == '-'
                        SEQS[i] = seq[1:j-1] * "." * seq[j:end]  
                        seq_len[i] += 1
                    end
                end
            end
            ncol = maximum(seq_len)  
        end
        is_insert = false  
        j += 1
    end
    SEQS
end

# Function to read A3M file and return IDS and SEQS
function _pre_reada3m(io)
    IDS, SEQS = MSA._pre_readfasta(io)  
    MSA._check_seq_and_id_number(IDS, SEQS) 
    try
        MSA._check_seq_len(IDS, SEQS) 
    catch
        SEQS = _add_inserts(SEQS) 
    end
    return IDS, SEQS 
end

# Example function to test the A3M parser
function test_a3m_parser()
    sequence_data = """
>example
ETESMKTVRIREKIKKFLGDRPRNTAEILEHINSTMRHGTTSQQLGNVLSKDKDIVKVGYIKRSGILSGGYDICEWATRNWVAEHCPEWTE
>1
----MRTTRLRQKIKKFLNERGeANTTEILEHVNSTMRHGTTPQQLGNVLSKDKDILKVATTKRGGALSGRYEICVWTLRP-----------
>2
----MDSQNLRDLIRNYLSERPRNTIEISAWLASQMDPNSCPEDVTNILEADESIVRIGTVRKSGMRLTDLPISEWASSSWVRRHE-----
>3
----MNSQNLRELIRNYLSERPRNTIEISTWLSSQIDPTNSPVDITSILEADDQIVRIGTVRKSGMRRSESPVSEWASNTWVKHHE-----
>4
--RDMDTEKVREIVRNYISERPRNTAEIAAWLNRH-DDGTGGSDVAAILESDGSFVRIGTVRTSGMTGNSPPLSEWATEKWIQHHER----
>5
-----RTRRLREAVLVFLEEKGnANTVEVFDYLNERFRWGATMNQVGNILAKDTRFAKVGHQ-RGQFRGSVYTVCVWALS------------
>6
-----RTKRLREAVRVYLAENGrSHTVDIFDHLNDRFSWGATMNQVGNILAKDNRFEKVGHVRD-FFRGARYTVCVWDLAS-----------
"""
    io = IOBuffer(sequence_data) 
    IDS, SEQS = _pre_reada3m(io) 

    println("IDS: ", IDS) 
    println("SEQS: ", SEQS)
end

# Run the test function
test_a3m_parser()

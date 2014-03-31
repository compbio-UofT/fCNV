#SAM bitewise flags
ALL_PROPERLY_ALIGNED = 0x2
IS_UNMAPPED = 0x4
MP_UNMAPPED = 0x8
#CIGAR operators
cigar_qr_op = 'HSIM=X'
cigar_db_op = 'NDM=X'
cigar_clip_op = 'HS'

def mapping_parser(m):
    '''
    Parse a read in SAM format, return a dictionary with filled in fields of interest.
    '''
    if isinstance(m, str):
        m = m.strip().split('\t')
    d = {}
    d['flag'] = int(m[1])   # flags
    d['chr'] = m[2]         # chr
    d['pos'] = int(m[3])    # pos
    d['mapq'] = int(m[4])   # mapping quality
    d['cigar'] = m[5]       # cigar string
    d['seq'] = m[9]         # sequence
    d['qual'] = m[10]       # sequencing quality
    
    return d

def parse_cigar_string(s):
    '''
    Parse given CIGAR string to a list of operators.
    '''
    res = []
    crt_len = ''
    i = 0
    while i < len(s):
        if str.isdigit(s[i]):
            crt_len += s[i]
        else:
            res.append([s[i], int(crt_len)])
            crt_len = ''
        i += 1  
    return res
    
def pile_up(read, data):
    '''
    Store allele evidence given by @read.
    '''
    #take only reads mapped in proper pair
    if read['flag'] & ALL_PROPERLY_ALIGNED == 0: return
    if read['mapq'] < 20: return
    
    pos_qr = 0
    pos_db = read['pos']
    op = parse_cigar_string(read['cigar'])
    for o in op:
        if o[0] == 'H': continue
        elif o[0] in 'SI': pos_qr += o[1]
        elif o[0] in 'ND': pos_db += o[1]
        elif o[0] in 'M=X':
            for i in range(o[1]):
                #if the read position is of sufficient quality, record this info
                if ord(read['qual'][pos_qr]) >= 33+20:
                    add_support(pos_db, read['seq'][pos_qr].upper(), data)
                pos_db += 1
                pos_qr += 1

def add_pos(pos, data):
    '''
    Add key @pos to @data.
    '''
    if not pos in data: data[pos] = dict()

def add_support(pos, nuc, data):
    '''
    Increase counter for nucleotide @nuc at position @pos, motify datastructure @data.
    Ignore positions @pos that are not in @data.
    '''
    if not pos in data: return #if we don't need info for this position
    if not nuc in 'ACGT': 
        print "Unrecognized symbol in SEQ:", nuc
        return
    try:
        data[pos][nuc] += 1
    except:
        data[pos].update({nuc:1})


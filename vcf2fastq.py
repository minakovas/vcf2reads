import os
import gzip
import random
import re
import argparse

GRCH38P14_CHR_PA = {'>NC_000001.11': 'chr1', '>NC_000002.12': 'chr2', '>NC_000003.12': 'chr3',
                    '>NC_000004.12': 'chr4', '>NC_000005.10': 'chr5', '>NC_000006.12': 'chr6',
                    '>NC_000007.14': 'chr7', '>NC_000008.11': 'chr8', '>NC_000009.12': 'chr9',
                    '>NC_000010.11': 'chr10', '>NC_000011.10': 'chr11', '>NC_000012.12': 'chr12',
                    '>NC_000013.11': 'chr13', '>NC_000014.9': 'chr14', '>NC_000015.10': 'chr15',
                    '>NC_000016.10': 'chr16', '>NC_000017.11': 'chr17', '>NC_000018.10': 'chr18',
                    '>NC_000019.10': 'chr19', '>NC_000020.11': 'chr20', '>NC_000021.9': 'chr21',
                    '>NC_000022.11': 'chr22', '>NC_000023.11': 'chrX', '>NC_000024.10': 'chrY'}


class VCF:
    """
    Creates object VCF

    Attributes:
    ----------
    metadata: list, contains strings, beginning from ##
    header: list with data column names
    samples: list with samples names(columns in data corresponding to the samples)
    data: list, each element contains data about one variation, elements correspond to the header
    """
    def __init__(self, path):
        _parse_vcf_out = self._parse_vcf(path)
        self.metadata = _parse_vcf_out[0]
        self.header = _parse_vcf_out[1]
        self.samples = _parse_vcf_out[2]
        self.data = _parse_vcf_out[3]

    @staticmethod
    def _parse_vcf(path):
        """
        Read a compressed vcf.vgz file

        Parameters
        ----------
        path: str, path to file.vcf.vgz

        Return
        -------
        tuple(metadata: list, header: list, samples: list, data: list)
        """
        with gzip.open(path, 'rb') as file:
            metadata = []
            line = file.readline()
            line = line.decode('utf-8')
            line = line.strip()
            while line.startswith('##'):
                metadata.append(line)
                line = file.readline()
                line = line.decode('utf-8')
                line = line.strip()
            header = line.split('\t')
            samples = header[header.index('FORMAT') + 1:]
            data = []
            for line in file:
                line = line.strip()
                line = line.decode('utf-8')
                data.append(line)
            data = [i.split('\t') for i in data]
            return metadata, header, samples, data


def parse_genome(path) -> tuple[dict, int]:
    """
    Read genome.fasta

    Parameters
    ----------
    path: str, path to genome.fasta

    Return
    ------
    genome: dict, keys are chromosomes/contigs, values are lists of strings from fasta with chromosome/contig sequence
    genome_line_len: int, number of nucleotides in one line in fasta file(constant length of line expected)

    Examples
    --------
    >>>parse_genome('genome.fasta')
    {'chr1':['NNNN....NNNN', 'NNNN....ATGC',...., 'TTGC....NNNN'], 'chr2': [...], ...}, 60
    """
    genome = {}
    with open(path, 'r') as file:
        seq_id = file.readline()
        seq_id = seq_id.strip()
        seq_id = GRCH38P14_CHR_PA[seq_id]
        seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                genome[seq_id] = seq
                seq_id = line
                seq_id = GRCH38P14_CHR_PA[seq_id]
                seq = []
            else:
                seq.append(line)
        genome_line_len = len(seq[0])
        genome[seq_id] = seq
    return genome, genome_line_len


def get_var_property(var_num: int, vcf: VCF, prop: str) -> str:
    """
    Get information from given column for given variation

    Parameters
    ----------
    var_num: int, index of variation in VCF.data(number of line in vcf file body)
    vcf: VCF, object storing information from vcf file
    prop: str, name of column in VCF.header

    Returns
    -------
    res: str, information from given column for given variation
    """
    prop_idx = vcf.header.index(prop)
    res = vcf.data[var_num][prop_idx]
    return res


def get_chr(var_num: int, vcf: VCF) -> str:
    """
    Returns chromosome for given variation

    Example
    -------
    >>>get_chr(10, vcf)
    'chrY'
    """
    return get_var_property(var_num, vcf, prop="#CHROM")


def get_pos(var_num: int, vcf: VCF) -> int:
    """
    Returns position in chromosome for given variation

    Example:
    >>>get_pos(10, vcf)
    1343635
    """
    pos = get_var_property(var_num, vcf, prop="POS")
    pos = int(pos)
    return pos


def get_ref(var_num: int, vcf: VCF) -> str:
    """
    Returns reference allele for given variation

    Example
    -------
    >>>get_ref(10, vcf)
    'A'
    """
    return get_var_property(var_num, vcf, prop="REF")


def get_alt(var_num: int, vcf: VCF) -> list:
    """
    Return list of alternative alleles for given variation

    Examples
    --------
    >>>get_alt(10, vcf)
    ['A']
    >>>get_alt(20, vcf)
    ['TTT', 'GCGC']
    """
    alt = get_var_property(var_num, vcf, prop="ALT")
    alt = alt.split(',')
    return alt


def get_format(var_num: int, vcf: VCF) -> list:
    """
    Return content of FORMAT column for given variation

    Example
    -------
    >>>get_format(10, vcf)
    ['GT', 'GC', 'DP', 'HQ']
    """
    frmt = get_var_property(var_num, vcf, prop="FORMAT")
    frmt = frmt.split(':')
    return frmt


def get_info(var_num: int, vcf: VCF) -> str:
    """
    Return content of INFO column for given variation

    Example
    -------
    >>>get_info(10, vcf)
    'NS=3;DP=14;AF=0.5;DB;H2'
    """
    return get_var_property(var_num, vcf, prop="INFO")


def get_sample_info(var_num: int, sample: str, vcf: VCF) -> list:
    """
    Return content of column, describing sample

    Example
    -------
    >>>get_sample_info(10, 'sample1', vcf)
    ['1|0', '1', '2', '0']
    """
    sample_info = get_var_property(var_num, vcf, prop=sample)
    sample_info = sample_info.split(':')
    return sample_info


def get_sample_gt(sample: str, var_num: int, vcf: VCF) -> str:
    """
    Returns genotype of the sample for given variation
    get_sample_info and get_format functions are used inside

    Example
    -------
    >>>get_sample_gt('sample1', 10, vcf)
    '0|1'
    """
    sample_info = get_sample_info(var_num, sample, vcf)
    sample_info_format = get_format(var_num, vcf)
    gt_idx = sample_info_format.index("GT")
    gt = sample_info[gt_idx]
    return gt


def get_called_variations(sample: str, vcf: VCF) -> list:
    """
    Return list of called variations for given sample
    get_sample_gt function is used inside
    Example
    -------
    >>>get_called_variations('sample1', vcf)
    [1, 2, 5, 7, ...]
    """
    sample_called_variants = []
    for i in range(len(vcf.data)):
        gt = get_sample_gt(sample, i, vcf)
        if gt == './.':  # not called
            continue
        else:
            sample_called_variants.append(i)
    return sample_called_variants


def get_variation_af(var_num: int, vcf: VCF) -> None or list:
    """
    Return list of alternative allele frequencies for given variation
    If AF not in INFO column in vcf, returns None
    If there is one alternative allele, returns [AF]
    If there are multiple alternative alleles, returns [AF1, AF2, ...]
    get_info function is used inside

    Example
    -------
    >>>get_variation_af(10, vcf) # one alternative allele
    [0.005]
    >>>get_variation_af(20, vcf) # three alternative alleles. for example G, GT, C
    [0.1, 0,5, 0.04]
    """
    info = get_info(var_num, vcf)
    pattern = r";AF=[^;]+"
    res = re.findall(pattern, info)
    if len(res) == 0:  # if not AF=x in INFO
        return None
    else:
        res = res[0].split("=")[1]
        res = res.split(",")
        af_lst = [float(i) for i in res]
    return af_lst


def isbetween(x, x_min, x_max) -> bool:
    """
    Returns True if x_min<x<x_max else False
    """
    return x_min <= x <= x_max


def is_variant_pass_af(gt: str, af_lst: list, af_min, af_max) -> bool:
    """
    Check if genotype pass by filtering by allele frequency
    If genotype is 0(reference allele) of 1, frequency of the first alternative allele is checked
    If genotype is 2, 3..., frequency of corresponding allele is checked
    isbetween function is used inside

    Parameters
    ----------
    gt: str, genotype, 0 - reference allele, 1 - alternative, if multiple alternative alleles, for example A, AT, G,
    1 - A, 2 - AT, 3 - G
    af_lst: list of allele frequencies for alternative alleles. Length of the list corresponds to the number of
    alternative alleles
    af_min: float, lower threshold for allele frequency
    af_max: float, higher threshold for allele frequency
    """
    if af_lst is None:  # if there was no AF=x in INFO column in vcf file
        return False
    if gt == '0' or gt == '1':
        return isbetween(af_lst[0], af_min, af_max)  # check first alternative allele
    else:  # gt = 2, 3...
        return isbetween(af_lst[int(gt)-1], af_min, af_max)  # check corresponding AF


def filter_variants_by_af(sample: str, sample_called_variants: list, vcf: VCF, af_min, af_max) -> list:
    """
    Filter list of variants by AF for given sample
    get_sample_gt, get_variation_af, is_variant_pass_af functions are used inside

    Parameters
    ----------
    sample: str, name of the sample to check
    sample_called_variants: list of variants to check, can be obtained using get_called_variations function
    vcf: VCF object, containing data from file.vcf
    af_min: float, lower threshold for allele frequency
    af_max: float, higher threshold for allele frequency
    """
    sample_called_variants_filtered = []
    for i in sample_called_variants:  # check each variation in sample_called_variants
        gt = get_sample_gt(sample, i, vcf)
        af_lst = get_variation_af(i, vcf)
        if is_variant_pass_af(gt, af_lst, af_min, af_max):
            sample_called_variants_filtered.append(i)
    return sample_called_variants_filtered


def find_position_in_chromosome(pos) -> tuple[int, int]:
    """
    Find line number and position in this line in the chromosome fasta file(i.e. in certain entry in fasta) for this
    position in the chromosome(from vcf). Numbering is 0-based.
    For example:
    ['NNNNNNNNN',
     'NNNANNNNN] for A line_number=0, pos_in_line=3

    !need global variable genome_line_len - number of nucleotides in one line in fasta

    Parameters
    ----------
    pos: int, position in chromosome(nucleotide number in linear sequence)

    Return
    ------
    line_number: int, number of line in fasta entry
    pos_in_line: int, number of nucleotide in this line
    """
    pos = int(pos) - 1  # from 1-based numbering in vcf to 0-based
    line_number = pos // genome_line_len
    pos_in_line = pos % genome_line_len
    return line_number, pos_in_line


def get_seq_chr(pos_start, pos_end, chr_lines: list) -> str:
    """
    Return chromosome region in interval [pos_start; pos_end)
    find_position_in_chromosome function is used inside

    Parameters
    ----------
    pos_start: int, start of the interval
    pos_end: int, end of the interval, not included!
    chr_lines: list, contains lines from fasta file
    """
    if pos_end < pos_start:
        raise ValueError("pos_end must be greater than pos_start")
    chr_start = find_position_in_chromosome(pos_start)
    chr_end = find_position_in_chromosome(pos_end)
    if chr_start[0] == chr_end[0]:
        seq = chr_lines[chr_start[0]][chr_start[1]:chr_end[1]]
        return seq
    else:
        seq = ''
        seq += chr_lines[chr_start[0]][chr_start[1]:]
        for line_num in range(chr_start[0]+1, chr_end[0]):
            seq += chr_lines[line_num]
        seq += chr_lines[chr_end[0]][:chr_end[1]]
        return seq


def replace_by_index(s: str, idx_left: int, old: str, new: str) -> str:
    """
    Replace old region to new in given sequence

    Parameters
    ----------
    s: str, some sequence
    idx_left: int, index of the left character of region to be replaced
    old: str, region to be replaced
    new: str, replacing region

    Examples
    --------
    >>>replace_by_index('aaTaa', 2, 'T', 'C')
    'aaCaa'
    >>>replace_by_index('aaTaa', 2, 'T', 'CCCCC')
    'aaCCCCCaa'
    >>>replace_by_index('aaTGTGTGaa', 2, 'TGTGTG', 'C')
    'aaCaa'
    """
    idx_right = idx_left + len(old)
    return s[:idx_left] + new + s[idx_right:]


def generate_random_subseq(seq, read_len) -> str:
    """
    selects a random fragment of length `read_len` from a given sequence `seq`
    For example generateRandomRead('ABCDEF', 3) --> 'ABC'|'BCD'|'CDE'|'DEF'
    """
    x = len(seq) - read_len
    start = random.randint(0, x)
    end = start + read_len
    read = seq[start:end]
    return read


def generate_read(pos, read_len, chr_lines: list, ref: str, alt=False) -> str:
    """
    Generate a random read involving the variation in the pos position
    Can generate insertions, deletions, SNPs and reads with reference allele
    If alt is given, ref is replaced by alt, if not, generate read with reference allele

    Algorithm:
    1) A fragment of the chromosome is taken from which the reed will be generated
       * new_pos - position of the first nucleotide of the reference allele in this fragment
       * left and right boundaries of this fragment are calculated
       * the fragment is obtained using get_seq_chr function
    2) if alt is given, reference is replaced by alt in the fragment using replace_by_index function
    3) Random read of given len is generated from the fragment using generateRandomRead function

    Parameters
    ----------
    pos: int, position of variation in chromosome
    read_len: int, length of read to generate
    chr_lines: list, containing lines from fasta of given chromosome
    ref: str, reference allele
    alt: str, alternative allele, replacing ref. Default False - generate read without replacing
    """
    new_pos = read_len-1
    left = pos - new_pos
    right = pos + (len(ref)-1) + new_pos
    seq = get_seq_chr(left, right+1, chr_lines)
    if alt:
        seq = replace_by_index(seq, new_pos, ref, alt)
    read = generate_random_subseq(seq, read_len)
    return read


def create_outdir(out_path):
    """
    Created directory for out_path if it is needed
    For example create_outdir('result/result1/test.fasta') creates directory ./result/result1/
    """
    directory = "/".join(out_path.split('/')[:-1]) # 'result/result1/test.fasta' --> 'result/result1
    if directory:
        if not os.path.exists(directory):
            os.makedirs(directory)


def generate_fasta_with_random_reads(vcf: VCF, sample, genome: dict, read_len, n, out_path, af_min, af_max) -> None:
    """
    Main function in the script
    Generate file.fasta in the given path with random
    create_outdir, get_called_variations, filter_variants_by_af, get_chr, get_pos, get_ref, get_sample_gt, get_alt,
    generate_read functions are used inside

    Parameters
    ----------
    vcf: VCF object
    sample: str, name of the sample
    genome: dict, containing chromosome(can be obtained by parse_genome function)
    read_len: int, length of reads to generate
    n: int, number of reads to generate
    out_path: str, path and name of file to write output fasta-file
    af_min: lower threshold for filtering by AF, if it is needed
    af_max: higher threshold for filtering by AF, if it is needed
    """
    create_outdir(out_path)
    vars_to_gen = get_called_variations(sample=sample, vcf=vcf)

    # filtering by AF, if af_min or af_max are not default
    if af_min != -1 or af_max != 2:
        vars_to_gen = filter_variants_by_af(sample=sample, sample_called_variants=vars_to_gen, vcf=vcf, af_min=af_min,
                                            af_max=af_max)
    if len(vars_to_gen) == 0:
        print(f"sample {sample} does not have called variations")
        return

    vars_to_gen = random.choices(vars_to_gen, k=n)
    file_out = open(out_path, 'w')
    # Generate read for each selected variation
    for i in vars_to_gen:
        chrom = get_chr(i, vcf)
        pos = get_pos(i, vcf)
        ref = get_ref(i, vcf)
        gt = get_sample_gt(sample, i, vcf)
        gt = int(gt)
        alt = get_alt(i, vcf)[gt - 1]
        if gt == 0:  # reference allele
            read = generate_read(pos=pos, read_len=read_len, chr_lines=genome[chrom], ref=ref)
        else:  # gt = 1, 2... alternative allele
            read = generate_read(pos=pos, read_len=read_len, chr_lines=genome[chrom], ref=ref, alt=alt)
        # write read to file
        read_id = f">Sample={sample}:Chr={chrom}:Pos={pos}:Ref={ref}:Alt={alt}:GT={gt}"
        print(read_id, file=file_out)
        print(read, file=file_out)
    file_out.close()


description = "This is a script for generating fasta file with random reads for some sample from reference genome " \
              "relying on VCF-file with optional filtering by allele frequency"

parser = argparse.ArgumentParser(description=description)

parser.add_argument("--genome", required=True, help="path to genome.fasta")
parser.add_argument("--vcf", required=True, help="path to file.vcf.vgz")
parser.add_argument("--read_len", required=True, type=int, help="length of read to generate")
parser.add_argument("--N", required=True, type=int, help="number of reads to generate")
parser.add_argument("--out", required=True, help="path to output file")
parser.add_argument("--sample", required=True, help="name of the sample in vcf")
parser.add_argument("--random_seed", default=42, help="random seed, default=42")
parser.add_argument("--af_min", default=-1, type=float, help="lower threshold for filtering by AF, if it is needed")
parser.add_argument("--af_max", default=2, type=float, help="higher threshold for filtering by AF, if it is needed")

args = parser.parse_args()

random.seed(args.random_seed)
genome, genome_line_len = parse_genome(args.genome)
vcf = VCF(args.vcf)

generate_fasta_with_random_reads(vcf=vcf,
                                 sample=args.sample,
                                 genome=genome,
                                 read_len=args.read_len,
                                 n=args.N,
                                 out_path=args.out,
                                 af_min=args.af_min,
                                 af_max=args.af_max)
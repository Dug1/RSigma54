import sys
import subprocess
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

class DisjointSets:
    def __init__(self, num):
        self.trees = [-1] * num

    def parent(self, num):
        if self.trees[num] < 0:
            return num
        else:
            parent = self.parent(self.trees[num])
            self.trees[num] = parent
            return parent 

    def join(self, a, b):
        parent_a = self.parent(a)
        parent_b = self.parent(b)
        if parent_a != parent_b:
            if self.trees[parent_a] > self.trees[parent_b]:
                self.trees[parent_b] += self.trees[parent_a]
                self.trees[parent_a] = parent_b
            else:
                self.trees[parent_a] += self.trees[parent_b]
                self.trees[parent_b] = parent_a

    def sets(self):
        sets = {}
        for i in range(len(self.trees)):
            index = self.parent(i)
            if index not in sets:
                sets[index] = []
            sets[index].append(i)
        return list(sets.values())

def meme(filenames, outdir):
    return subprocess.run(["meme-meme"] + filenames + ["-dna","-oc", outdir])

def rpstblastn(input_file_name, output_file_name, database_path):
    return subprocess.run(["rpstblastn", "-query", input_file_name, "-out", output_file_name, "-db", database_path, "-outfmt", "5"])

def blastp(queryfile, db, outfile, e=10):
    return subprocess.run(["blastp", "-query", queryfile, "-db", db, "-out", outfile, "-evalue", str(e), "-outfmt", "5"])

def makeblastdb(f, db_type):
    return subprocess.run(["makeblastdb", "-in", f, "-dbtype", db_type])

def extract_locus(string):
    return string.split(" ")[0].split("|")[1]

def extract_id(string):
    return string.split(" ")[0].split("|")[0]

def extract_full_id(string):
    return string.split(" ")[0]

def filter_matches(filename, e):
    """ Parses the the matches results of rpstblastn and save every match 
    above a certain evalue"""
    transcription_factors = []
    with open(filename) as f:
        try:
            for record in NCBIXML.parse(f):
                best = e
                if record.alignments:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect < best:
                                best = hsp.expect

                if best < e:
                    transcription_factors.append(extract_locus(record.query))
        except ValueError as e:
            return None

    return transcription_factors

def find_matches(filename, e):
    """ Find the best match for every query protein above a certain e-value, seperated by ids
    above a certain evalue"""
    best_matches = {}
    with open(filename) as f:
        try:
            for record in NCBIXML.parse(f):
                best = {}
                if record.alignments:
                    for alignment in record.alignments:
                        genome = extract_id(alignment.hit_def)
                        locus = extract_locus(alignment.hit_def)
                        
                        best_value = e
                        for hsp in alignment.hsps:
                            if hsp.expect < best_value:
                                best_value = hsp.expect
                        
                        if genome not in best:
                            best[genome] = []
                        
                        best[genome].add((locus, best_value))

                    best_matches[extract_full_id(record.query)] = best

        except ValueError as e:
            return None

    return best_matches 

def extract_genes(filename, fmt, outname, outfmt, locuses=None):
    """ Extra all proteins coded in the file """

    genome = SeqIO.parse(filename, fmt)
    proteins = []
    
    if fmt == "genbank":
        for record in genome:
            for feature in filter(lambda f: f.type == "CDS", record.features):
                qualifiers = feature.qualifiers
                if locuses is None or qualifiers["locus_tag"][0] in locuses:
                    name = qualifiers["gene"][0] if "gene" in qualifiers else ""
                    desc = qualifiers["product"][0] if "product" in  qualifiers else ""
                    protein = SeqRecord(Seq(qualifiers["translation"][0], IUPAC.protein), 
                        id="{0}|{1}".format(record.id, qualifiers["locus_tag"][0]),
                        name=name,
                        dbxrefs=qualifiers["db_xref"],
                        description=desc)

                    proteins.append(protein)
    elif fmt == "fasta":
        if locuses == None:
            proteins = genome
        else:
            for seq in genome:
                if extract_locus(seq.id) in locuses:
                    proteins.append(seq)

    SeqIO.write(proteins, outname, outfmt)

def build_file_lookup(files):
    file_map = {}
    for fname in files:
        for seq in SeqIO.parse(fname, "fasta"):
            file_map[extract_id(seq.id)] = fname

    return file_map

def merge_entries(job_id, requests):
    connections = []
    locus = set()
    for request in requests:
        locuses = ["{0}|{1}".format(obj["file"], obj["locuses"][0]) for obj in request]
        for i in range(len(locuses)):
            locus.add(locuses[i])
            for j in range(1, i + 1):
                connections.append((locuses[i], locuses[j]))
    
    locus = list(locus)
    pairs = zip(locus, range(len(locus)))
    lookup = dict(pairs)
    ds = DisjointSets(len(locus))
    for first, second in connections:
        ds.join(lookup[first], lookup[second])
        
    sets = ds.sets()
    print(sets)
    for i in range(len(sets)):
        data = []
        for index in sets[i]:
            iden = locus[index].split("|")
            data.append({"file":iden[0], "locuses":[iden[1]]})
            result_name = "{0}-{1}-memereq.json".format(job_id, i)
            with open(result_name, "w") as f:
                f.write(json.dumps(request))

def find_transcription_factors(query_id, filename, output_file, profiledb, e):
    """ Finds all transcription factors which match the [profiledb] 
    and saves them to a fasta file"""

    temp_file = "{0}-tf.fsa".format(query_id)
    extract_genes(filename + ".gbk", "genbank", temp_file, "fasta")
    outfile = "{0}-tf.xml".format(query_id)
    print(rpstblastn(temp_file, outfile, profiledb))
    locuses = filter_matches(outfile, e)
    print(locuses)
    extract_genes(filename + ".gbk", "genbank", output_file + ".fsa", "fasta", locuses)

def generate_match_matrix(job_id, original, test_against):
    file_lookup = build_file_lookup(test_against)
    genes = []
    for filename in test_against:
        genes += SeqIO.parse(filename, "fasta")

    filename = "{0}-congregate.fsa".format(job_id)
    SeqIO.write(genes, filename, "fasta")
   
    makeblastdb(filename, "prot")
    out_file = "{0}-ortho".format(job_id)
    blastp(original, filename, out_file, 1)
    best_locuses = find_best_matches(out_file, 1)

    print(len(best_locuses))

    condensed_matrix = dict()
    for orig, data in best_locuses.items():
        file_map = dict()
        for iden, tup in data.items():
            f = file_lookup[iden]
            locus, e = tup
            if f not in file_map or file_map[f][1] > e:
                file_map[f] = tup

        condensed_matrix[orig] = file_map

    return condensed_matrix

def find_orthologs_designated(job_id, query, original_file, test_against):
    """Finds orthologs by buiding a lookup matrix between original_file and test against and
    verifying the result by comparing matches on all other genes."""
    matrix = generate_match_matrix(job_id, original_file, test_against)
    for element in SeqIO.parse(query, "fasta"):
        
    print(matrix)

def find_orthologs_old(job_id, test_against_matrix):
    """ Finds the orthologs of in a given set of each query set of genes against a list
    of others. Uses multiple queries of the blastp"""

    requests = []
    for query, data in test_against_matrix.items():
        print("\n Processing matches for {0}".format(query))

        original_file, test_against = data
        makeblastdb(original_file, "prot")
        orthologs = dict()
        file_lookup = build_file_lookup(test_against)

        genes = []
        for filename in test_against:
            genes += SeqIO.parse(filename, "fasta")

        filename = "{0}-congregate.fsa".format(job_id)
        SeqIO.write(genes, filename, "fasta")

        makeblastdb(filename, "prot")
        out_file = "{0}-ortho".format(job_id)
        blastp(query, filename, out_file)
        best_locuses = find_best_matches(out_file, 1)

        if best_locuses == None:
            continue

        locus_list = []
        for locus, matches in best_locuses.items():
            for genome, iden in  matches.items():
                locus_list.append(iden[0])

        query_file = "{0}-query.fsa".format(job_id)
        extract_genes(filename, "fasta", query_file, "fasta",  set(locus_list))

        blastp(query_file, original_file, out_file)
        reverse_matrix = find_best_matches(out_file, 1)
        if reverse_matrix == None:
            continue

        for target_locus, matches in best_locuses.items():
            for genome, loc_pair in matches.items():
                tag = "{0}|{1}".format(genome, loc_pair[0])
                ortholog = list(reverse_matrix[tag].items())[0]
                other_locus = "{0}|{1}".format(ortholog[0], ortholog[1][0])
                if target_locus == other_locus:
                    if not target_locus in orthologs:
                        orthologs[target_locus] = []
                    orthologs[target_locus].append(tag)
    
        for original_locus, orthologs in orthologs.items():
            request = [{"file":original_file, "locuses":[extract_locus(original_locus)]}]
            for ortho_id in orthologs:
                gen_id, locus = tuple(ortho_id.split("|"))
                request.append({"file": file_lookup[gen_id], "locuses":[locus]})
            
            requests.append(request)

    return requests

def extract_upstream(indicies, genome, amount, overlap, min_length=8):
    """ For each index in [indicies], extract [amount] codons before it or
        up until the previous genes is overlap is false """

    records = []
    prev_end = -1
    index = 0
    for feature in filter(lambda f: f.type == "CDS", genome.features):
        if index in indicies:
            end = int(feature.location.start)
            start = max(end - amount, 0)
            if not overlap:
                start = max(start, prev_end)

            if (end - start) > min_length:
                upstream = genome[start:end]
                upstream.id =  "{0}|{1}".format(genome.id, feature.qualifiers["locus_tag"][0])
                records.append(upstream)

        index += 1
        prev_end = int(feature.location.end)

    return records


def extract_upstream_for_meme(genomes, locuses, upstream, radius, overlap):
    """ Given genomes and locuses, returns all [upstream] number
        of codons before the gene and its [radius] number of
        surrounding genes. [overlap] determines whether or not to allow the upstream to
        bleed into the previous CDS."""

    records = []
    for genome in genomes:
        feature_len = len(genome.features)

        index = 0
        locations = set()
        for feature in filter(lambda f: f.type == "CDS", genome.features):
            locus = feature.qualifiers["locus_tag"][0] 
            if locus in locuses:
                locations.add(index)
                for i in range(index - radius, index + radius):
                    locations.add(i)
            
            index += 1

        print(locations)
        records += extract_upstream(locations, genome, upstream, overlap)

    return records

    
def find_tfbs(job_id, filename, outdir, upstream, radius, overlap):
    """ Finds the the best motif using using the query specified in [filename]
        and outputs the files into the directory specified in [outdir]"""

    records = []
    with open(filename) as f:
        jobs = []
        for req in json.loads(f.read()):
            print(req["file"])
            genomes = list(SeqIO.parse(req["file"], "genbank"))
            records += extract_upstream_for_meme(genomes, req["locuses"], upstream, radius, overlap)

    temp_file = "{0}-meme.fsa".format(job_id)
    SeqIO.write(records, temp_file, "fasta")
    meme([temp_file], outdir)
        
        
#find_transcription_factors("1", "511145", "results", "domain-matching/db/Sigma.pin", 1)
#find_tfbs(1, "meme/request.json", "output", 150, 0, False)
#jobs = frozenset(["211586", "225849", "318167", "319224", "323850", "325240", "326297", "351745", "392500", "398579", "425104", "458817", "60480", "60481", "94122"])
#matrix = dict()
#for job in jobs:
#    test_against = ["meme/{0}.fsa".format(other) for other in jobs - set([job])]
#    matrix["query/" + job + ".fsa"] = ("meme/{0}.fsa".format(job), test_against)
#merge_entries(1, find_orthologs(1, matrix))
designated = "meme/225849.fsa"
query = "query/225849.fsa"
jobs = frozenset(["211586", "318167", "319224", "323850", "325240", "326297", "351745", "392500", "398579", "425104", "458817", "60480", "60481", "94122"])
test_against = ["meme/{0}.fsa".format(other) for other in jobs]
find_orthologs_designated(1, query, designated, test_against)

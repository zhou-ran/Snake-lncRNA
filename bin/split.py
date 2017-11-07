from Bio import SeqIO
import sys
import os


Sequence , Num  ,Path = sys.argv[1:]
filetype = Sequence.split('.')[-1]
#os.mkdir(Path)
#os.chdir(Path)
#print Sequence
""""
usage:

python split.py yourfile num

your file's suffix must be 'fasta' or 'fastq'
"""
#print Sequence , Num , filetype

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


record_iter = SeqIO.parse(open(Sequence),filetype)
for i, batch in enumerate(batch_iterator(record_iter, int(Num))):
    dirname = Path + '/split-%i' % (i + 1)
    if os.path.isdir(dirname):
        pass
    else:
        os.makedirs(dirname)
    filename = 'split.' + filetype
    handle = open(dirname + "/" + filename, "w")
    count = SeqIO.write(batch, handle, filetype)
    handle.close()
    print("Wrote %i records to %s" % (count, filename))
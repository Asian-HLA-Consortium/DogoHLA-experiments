"""Frozen strict unphased diploid endpoints for named four-field and core identities."""
import re

def full_label(x):
 x=x.removeprefix('HLA-')
 return x if re.fullmatch(r'[A-Z0-9]+\*\d+:\d+:\d+:\d+[A-Z]*',x) else None

def decode_slots(line,aliases):
 v=line.rstrip().split('\t');copies=int(v[1]);tokens=[v[2] if float(v[4])>0 else '.'] if copies in [1,2] else ['.']
 tokens=[tokens[0],tokens[0] if copies==1 else v[5] if copies==2 and float(v[7])>0 else '.']
 out=[]
 for token in tokens:
  labels=set();cores=set();label_unknown=False;core_unknown=False
  for name in token.split(','):
   a=aliases.get(name,aliases.get(name.removeprefix('HLA-')))
   if a is not None:
    ls={full_label(x) for x in a['genomic_labels']};label_unknown|=not ls or None in ls;labels.update(x for x in ls if x)
    cs=set(a.get('core_ids',[]));core_unknown|=not cs;cores.update(cs)
   else:
    label=full_label(name);label_unknown|=label is None
    if label:labels.add(label)
    # A named allele alone does not establish an arbitrary assembly-annotation core.
    core_unknown=True
  out.append(dict(labels=sorted(labels),cores=sorted(cores),label_unknown=label_unknown,core_unknown=core_unknown))
 return out

def strict_pair(predicted,truth,kind):
 key='labels' if kind=='four_field' else 'cores';unknown='label_unknown' if kind=='four_field' else 'core_unknown'
 # No ambiguity credit, no exclusion of incorrect/no-call predictions from denominator.
 eligible=len(truth)==2 and all(len(x)==1 for x in truth)
 called=len(predicted)==2 and all(len(x[key])==1 and not x[unknown] for x in predicted)
 return dict(eligible=eligible,called=called,correct=eligible and called and sorted(x[key][0] for x in predicted)==sorted(x[0] for x in truth))

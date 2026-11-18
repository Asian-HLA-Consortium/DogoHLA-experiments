"""Shared numeric four-field and exact provisional-core endpoints; no padding."""
import itertools,re
PATTERN=re.compile(r'(?:HLA-)?([A-Z0-9]+\*\d+:\d+:\d+:\d+)([NLSCAQ]?)')
def numeric_four(label):
 m=PATTERN.fullmatch(label)
 return m[1] if m else None

def truth_slots(copies):
 out=[]
 for copy in copies:
  names=copy['genomic_labels'];values={numeric_four(n) for n in names}
  # A shorter alternative is unresolved, rather than silently discarded.
  out.append({'labels':sorted(x for x in values if x),'unknown':not names or None in values})
 return out

def named_pair(prediction,copies):
 truth=truth_slots(copies)
 eligible=len(truth)==2 and all(not x['unknown'] and len(x['labels'])==1 for x in truth)
 slots=[]
 for p in prediction:
  values={numeric_four(x) for x in p['labels']}
  slots.append({'labels':sorted(x for x in values if x),'unknown':p['label_unknown'] or not values or None in values})
 called=len(slots)==2 and all(not x['unknown'] and len(x['labels'])==1 for x in slots)
 correct=eligible and called and sorted(x['labels'][0] for x in slots)==sorted(x['labels'][0] for x in truth)
 return {'eligible':eligible,'called':called,'correct':correct,'truth_expression_suffix_present':any(PATTERN.fullmatch(n) and PATTERN.fullmatch(n)[2] for c in copies for n in c['genomic_labels'])}

def unmatched_cores(prediction,copies):
 expected=[c['core_id'] for c in copies if not c['genomic_labels']]
 if len(expected)>2 or len(prediction)!=2:raise ValueError('Expected diploid records')
 correct=max((sum(not prediction[j]['core_unknown'] and prediction[j]['cores']==[core] for j,core in zip(order,expected)) for order in itertools.permutations(range(2),len(expected))),default=0)
 return {'eligible_copies':len(expected),'exact_copies_recovered':correct}

def canonical_slots(call):
 # Canonical Dogo/FuFi records can be converted to the same named endpoint.
 out=[]
 for slot in call.get('resolutions',{}).get('4',[]):
  out.append({'labels':slot.get('alleles',[]),'label_unknown':slot.get('unresolved',True) or slot.get('harmonization_unresolved',False),'cores':[],'core_unknown':True})
 return out

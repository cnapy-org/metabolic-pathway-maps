small = ['adp', 'atp']

def findSmallSpeciesAssociatedWithReaction(cmod, small):
    out = {}
    for r_ in cmod.reactions:
        i = set(r_.getSpeciesIds()).intersection(set(small))
        if len(i) > 0:
            out[r_.getId()] = list(i)
            
    return out
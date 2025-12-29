from collections import defaultdict

class IDGenerator: 
    """Generate unique graph IDs with semantic prefixes.

    The generator maintains separate counters per prefix and returns IDs
    such as ``m0001`` for metabolites, ``e0001`` for enzymes, or ``p0001``
    for pathways. Prefixes are inferred from the requested kind.

    Attributes
    ----------
    counters : collections.defaultdict
        Mapping from single-letter prefix to the next integer counter
        value used to build IDs.
    """

    def __init__(self):
        """Initialize the ID generator with empty prefix counters."""
        self.counters = defaultdict(int)        
 

    def _prefix(self, kind: str) -> str:  
        """Map an entity kind string to a single-letter ID prefix.

        Parameters
        ----------
        kind : str
            Kind of object to identify (e.g. "metabolite", "enzyme",
            "pathway", "interaction").

        Returns
        -------
        str
            One-letter prefix code: ``"m"`` for metabolites, ``"e"``
            for enzymes or gene products, ``"p"`` for pathways,
            ``"a"`` for anchors, ``"i"`` for interactions, or
            ``"n"`` for unrecognized kinds.
        """       
        k = (kind or "").lower()
        if k in {"enzyme", "geneproduct", "gene", "protein", "e"}: return "e"
        if k in {"metabolite", "compound", "smallmolecule", "m"}:  return "m"
        if k in {"pathway"}:                                       return "p"
        if k in {"anchor", "a"}:                                   return "a"
        if k in {"interaction", "edge", "i"}:                      return "i"
        return "n"

    def new(self, kind: str) -> str:
        """Create a new unique ID for a given kind of object.

        Parameters
        ----------
        kind : str
            Kind of object to identify; passed to ``_prefix`` to
            determine the prefix.

        Returns
        -------
        str
            Newly generated identifier in the form ``<prefix><NNNN>``,
            where the numeric part is zero-padded to four digits and
            monotonically increases per prefix.
        """
        p = self._prefix(kind)
        self.counters[p] += 1
        return f"{p}{self.counters[p]:04d}"  # p=prefix + 4 digits with leading 0's 


# initiation of an instance for this class
idgenerator = IDGenerator()
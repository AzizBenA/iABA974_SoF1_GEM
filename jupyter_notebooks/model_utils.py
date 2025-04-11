
from collections import defaultdict, Counter
from collections.abc import Sequence
import re
from typing import TypeVar, Literal
import warnings

import mergem
from cobra.core import Model, Reaction, Metabolite, Group, Gene
from cobra.core.formula import Formula
from pandas import Series, DataFrame, concat
from tabulate import tabulate


def prep_sof_model_for_merge ( orig_sof:Model ) -> Model :
    sof = orig_sof.copy()
    # Sof1 model uses non-standard "C_c" for compartment.
    for met in sof.metabolites :
        match met.compartment :
            case "c"|"C_c" : met.compartment = "c"
            case "C_p" : met.compartment = "p"
            case "C_e" : met.compartment = "e"
            case _ : Exception(f"Cannot handle compartment {met.compartment}")
    return sof


def prep_xau_model_for_merge ( orig_xau:Model ) -> Model :
    # Xau model uses IDs with compartments that are unrecognizable by mergem.
    xau = orig_xau.copy()
    for met in xau.metabolites :
        m = re.search(r"^([CG][0-9]+)", met.id)
        if m is None :
            raise Exception(f"Cannot deal with metabolic ID = {met.id}")
        met.id = m.group(1)
    return xau


def compare_formula_similarity ( formula_1: str, formula_2: str, strict: bool = False ) -> bool :
    """Compare two formulas and apply some similarity rules to elements."""
    formula_1_elements = Formula(formula_1).elements
    formula_2_elements = Formula(formula_2).elements
    if not strict :
        if "H" in formula_1_elements: del formula_1_elements["H"]
        if "H" in formula_2_elements: del formula_2_elements["H"]
    return formula_1_elements == formula_2_elements


# TODO: There might be a compare reaction equation similarity as well. Maybe it can ignore co-factors.


def collect_kegg_reactions_from_model ( model:Model, include_similar:bool = False ) -> set[str] :
    """
    Collect a set of all Kegg reactions a model has, either from the ID or from a kegg annotation.
    @param include_similar
        If true, will use mergem to include Kegg IDs from reactions that are similar to those
        inside the model, even if those similar Kegg IDs to not exist in the model.
    """
    reactions = list()
    for reac in model.reactions :
        poss_reacs = set()  # possible reactions that may be added
        poss_reacs.add(reac.id)
        for annot_key in ["kegg","kegg.reaction"] :
            if annot_key in reac.annotation :
                if isinstance(reac.annotation[annot_key],list) :
                    poss_reacs.update(reac.annotation[annot_key])
                elif isinstance(reac.annotation[annot_key],str) :
                    poss_reacs.add(reac.annotation[annot_key])
        filt_reacs = { m.group(1) for id in poss_reacs if (m := re.match(r"(R[0-9]+)", id)) is not None }
        reactions.extend(filt_reacs)

    # inside the same model reaction, we allow duplicate Kegg ID, like when name and
    # annotation are equal, but we don't want multiple different model reactions
    # to use the same Kegg IDs.
    dup_reactions = [ f"{k}_x{v}" for k,v in Counter(reactions).items() if v > 2 ]
    if len(dup_reactions) > 0 :
        warnings.warn(f"There are {len(dup_reactions)} duplicate reactions in model {model.name}.")
    reaction_set = set(reactions)  # From this point we don't care about duplicates anymore.

    if include_similar is True :
        for rid in reaction_set :
            similar_rids = get_similar_reactions_from_mergem(rid)
            reaction_set.update(similar_rids)

    return reaction_set


def collect_kegg_metabolites_from_model ( model:Model, include_similar:bool = False ) -> set[str] :
    """
    Collect a set of all Kegg compounds a model has, either from the ID or from a kegg annotation.
    @param include_similar
        If true, will use mergem to include Kegg IDs from metabolites that are similar to those
        inside the model, even if those similar Kegg IDs to not exist in the model.
    """
    metabolites = list()
    for met in model.metabolites :
        poss_mets = set()  # possible metabolites that may be added
        poss_mets.add(met.id)
        for annot_key in ["kegg","kegg.compound","kegg.drug","kegg.glycan"] :
            if annot_key in met.annotation :
                if isinstance(met.annotation[annot_key],list) :
                    poss_mets.update(met.annotation[annot_key])
                elif isinstance(met.annotation[annot_key],str) :
                    poss_mets.add(met.annotation[annot_key])
        filt_mets = { m.group(1) for id in poss_mets if (m := re.match(r"(C[0-9]+)", id)) is not None }
        metabolites.extend(filt_mets)

    # inside the same model metabolite, we allow duplicate Kegg ID, like when name and
    # annotation are equal, but we don't want multiple different model metabolites
    # to use the same Kegg IDs.
    dup_metabolites = [ f"{k}_x{v}" for k,v in Counter(metabolites).items() if v > 2 ]
    if len(dup_metabolites) > 0 :
        warnings.warn(f"There are {len(dup_metabolites)} duplicate metabolites in model {model.name}.")
    metabolite_set = set(metabolites)

    if include_similar is True :
        for mid in metabolite_set :
            similar_mids = get_similar_metabolites_from_mergem(mid)
            metabolite_set.update(similar_mids)

    return metabolite_set


def get_similar_reactions_from_mergem ( reaction_id:str, source_db:str|None = None ) -> set[str] :
    """
    Checks the mergem database to find all similar reactions (including this one)
    that are sourced from the same database.
    If `source_db` is specified, it will get the entries from that DB instead of using
    the same DB as the input.
    Example:
        ("R01067") -> {"R01067","R01830"}
    TODO: Method is too shallow, consider removing it.
    """
    return _get_similar_entries_from_mergem( "REACTION", reaction_id, source_db )


def get_similar_metabolites_from_mergem ( metabolite_id:str, source_db:str|None = None ) -> set[str] :
    """
    Checks the mergem database to find all similar metabolites (including this one)
    that are sourced from the same database.
    If `source_db` is specified, it will get the entries from that DB instead of using
    the same DB as the input.
    Example:
        ("C00085") -> {"C00085","C05345"}
        ("f6p") -> {"f6p","f6p_B"}
    TODO: Method is too shallow, consider removing it.
    """
    return _get_similar_entries_from_mergem( "METABOLITE", metabolite_id, source_db )


def _get_similar_entries_from_mergem (
    entry_type:Literal["METABOLITE","REACTION"], input_id:str, source_db:str|None = None
) -> set[str] :
    prop: dict
    match entry_type :
        case "METABOLITE" :
            univ_id = mergem.map_metabolite_univ_id(input_id)
            met_prop = mergem.get_metabolite_properties(univ_id)
            if not isinstance(met_prop,dict) :
                return set()
            prop = met_prop
        case "REACTION" :
            univ_id = mergem.map_reaction_univ_id(input_id)
            reac_prop = mergem.get_reaction_properties(univ_id)
            if not isinstance(reac_prop,dict) :
                return set()
            prop = reac_prop

#     # Find from which database the ínput ID is sourced from.
#     source_db: str|None = None
#     mid: str
#     for db_mid in met_prop["ids"] :
#         (db, mid) = db_mid.split( sep=":", maxsplit=1 )
#         if db == metabolite_id :
#             source_db = db
#             break
#     if source_db is None :
#         raise Exception("Internal error. Source database from metabolite cannot be found.")

    # Find from which database the ínput ID is sourced from if no fixed DB is provided.
    if source_db is None :
        source_dbs:list[str] = [
            ss[0]
            for db_mid in prop["ids"]
            if (ss := db_mid.split( sep=":", maxsplit=1 )) and ss[1] == input_id
        ]
        if len(source_dbs) != 1 :
            print(f"Internal error. Mergem mapped {input_id} to an entry, but that entry does not include it.")  # TODO use logger
            return set()
        source_db = source_dbs[0]

    # Collect all IDs that use the same source database.
    similar_ids:list[str] = [
        ss[1]
        for db_mid in prop["ids"]
        if (ss := db_mid.split( sep=":", maxsplit=1 )) and ss[0] == source_db
    ]
    return set(similar_ids)


CobraEntry = TypeVar('CobraEntry', Reaction, Metabolite, Group, Gene)

def search_index_by_annotation_value ( entries:Sequence[CobraEntry], value:str ) -> list[int] :
    """
    Given a sequence of cobra `entries`, typically `cobra_model.reactions` or `cobra_model.metabolites`,
    return the indexes that match the `value` in some form (id or annotation).
    """
    found_entries = list()
    for i,entry in enumerate(entries) :
        found_annot_vals = set()
        found_annot_vals.add(entry.id)  # ID is also an annotation
        for annot_key,annot_val in entry.annotation.items() :
            if isinstance(annot_val,list) :
                found_annot_vals.update(annot_val)
            elif isinstance(annot_val,str) :
                found_annot_vals.add(annot_val)
        if value in found_annot_vals :
            found_entries.append(i)
    return found_entries


class ModelSummary :

    rows: list[dict]

    def __init__ ( self ) -> None :
        self.rows = list()

    def add( self, model:Model, alt_name: str|None = None ) -> "ModelSummary" :
        self.rows.append({
            "Model": model.name if alt_name is None else alt_name,
            "Metabolites": len(model.metabolites),
            "Reactions": len(model.reactions),
            "Genes": len(model.genes),
            "Groups": len(model.groups),
        })
        return self

    def output (self) -> str :
        return tabulate( tabular_data=self.rows, headers="keys" )


class ModelAnnotationStats :

    serieses: list[Series] = list()

    def __init__ ( self ) -> None :
        self.serieses = list()

    def add( self, model:Model, alt_name: str|None = None ) -> "ModelAnnotationStats" :
        stats = Counter()
        for entry_key in ("genes","groups","metabolites","reactions") :
            entries = getattr(model, entry_key)
            for entry in entries :
                stats[f"{entry_key}.id.{type(entry.id).__name__}"] += 1
                stats[f"{entry_key}.name.{type(entry.name).__name__}"] += 1
                for a_key, a_val in entry.annotation.items() :
                    stats[f"{entry_key}.a.{a_key}.{type(a_val).__name__}"] += 1
                for n_key, n_val in entry.notes.items() :
                    stats[f"{entry_key}.n.{n_key}.{type(n_val).__name__}"] += 1
        self.serieses.append(
            Series( data=stats, name=(model.name if alt_name is None else alt_name) )
        )
        return self

    def as_dataframe (self) -> DataFrame :
        return concat( self.serieses, axis="columns" ).sort_index()

    def output (self) -> str :
        return self.as_dataframe().to_string()


import os
import pandas as pd
import numpy as np

# Locate the database file relative to the module's path
hugo_db_path = os.path.join(os.path.dirname(__file__), "hugo_export.txt")
hugo_db = pd.read_csv(hugo_db_path, sep="\t", low_memory=False)


def lookup_hugo_gene_symbol(gene_symbol: str) -> dict:
    ret_cols = [
        "hgnc_id",
        "ensembl_gene_id",
        "uniprot_ids",
        "symbol",
        "alias_symbol",
        "prev_symbol",
    ]
    matches = hugo_db.loc[hugo_db.symbol == gene_symbol, ret_cols]
    if matches.shape[0] == 0:
        mask = hugo_db.prev_symbol.str.contains(
            rf"\b{gene_symbol}\b", regex=True
        )
        mask[mask.isna().values] = False
        matches = hugo_db.loc[mask, ret_cols]
        if matches.shape[0] == 0:
            mask = hugo_db.alias_symbol.str.contains(
                rf"\b{gene_symbol}\b", regex=True
            )
            mask[mask.isna().values] = False
            matches = hugo_db.loc[mask, ret_cols]
            if matches.shape[0] == 0:
                raise ValueError(
                    "Gene symbol does not match anything in database."
                )
    matches = np.array(matches)[~matches.isna().values]
    out_dic = {
        "hugo_id": matches[0],
        "ensembl_id": matches[1],
        "uniprot_id": matches[2],
        "gene_symbols": [
            item
            for sublist in [m.split("|") for m in matches[3:]]
            for item in sublist
        ],
    }
    return out_dic


def lookup_hugo_gene_id(gene_id: int) -> dict:
    gene_id = f"HGNC:{gene_id}"
    ret_cols = [
        "hgnc_id",
        "ensembl_gene_id",
        "uniprot_ids",
        "symbol",
        "alias_symbol",
        "prev_symbol",
    ]
    matches = hugo_db.loc[hugo_db.hgnc_id == gene_id, ret_cols]
    if matches.shape[0] == 0:
        raise ValueError("Gene id does not match anything in database.")
    matches = np.array(matches)[~matches.isna().values]
    out_dic = {
        "hugo_id": matches[0],
        "ensembl_id": matches[1],
        "uniprot_id": matches[2],
        "gene_symbols": [
            item
            for sublist in [m.split("|") for m in matches[3:]]
            for item in sublist
        ],
    }
    return out_dic


def lookup_uniprot_accession(accession: str) -> dict:
    ret_cols = [
        "hgnc_id",
        "ensembl_gene_id",
        "uniprot_ids",
        "symbol",
        "alias_symbol",
        "prev_symbol",
    ]
    mask = hugo_db.uniprot_ids.str.contains(rf"\b{accession}\b", regex=True)
    mask[mask.isna().values] = False
    matches = hugo_db.loc[mask, ret_cols]
    if matches.shape[0] == 0:
        raise ValueError("Uniprot does not match anything in database.")
    matches = np.array(matches)[~matches.isna().values]
    out_dic = {
        "hugo_id": matches[0],
        "ensembl_id": matches[1],
        "uniprot_id": matches[2],
        "gene_symbols": [
            item
            for sublist in [m.split("|") for m in matches[3:]]
            for item in sublist
        ],
    }
    return out_dic


def lookup_ensembl_gene_symbol(gene_symbol: str) -> dict:
    ret_cols = [
        "hgnc_id",
        "ensembl_gene_id",
        "uniprot_ids",
        "symbol",
        "alias_symbol",
        "prev_symbol",
    ]
    matches = hugo_db.loc[hugo_db.ensembl_gene_id == gene_symbol, ret_cols]
    if matches.shape[0] == 0:
        mask = hugo_db.prev_symbol.str.contains(
            rf"\b{gene_symbol}\b", regex=True
        )
        mask[mask.isna().values] = False
        matches = hugo_db.loc[mask, ret_cols]
        if matches.shape[0] == 0:
            mask = hugo_db.alias_symbol.str.contains(
                rf"\b{gene_symbol}\b", regex=True
            )
            mask[mask.isna().values] = False
            matches = hugo_db.loc[mask, ret_cols]
            if matches.shape[0] == 0:
                raise ValueError(
                    "Gene symbol does not match anything in database."
                )
    matches = np.array(matches)[~matches.isna().values]
    out_dic = {
        "hugo_id": matches[0],
        "ensembl_id": matches[1],
        "uniprot_id": matches[2],
        "gene_symbols": [
            item
            for sublist in [m.split("|") for m in matches[3:]]
            for item in sublist
        ],
    }
    return out_dic

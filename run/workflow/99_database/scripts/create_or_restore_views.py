import sqlite3
import argparse

def create_views(con):
    cur = con.cursor()
    ### View for prediction scores - v_predictions
    cur.execute('DROP VIEW IF EXISTS v_predictions')
    
    scores_distinct = cur.execute(r'''
        SELECT DISTINCT prediction_tool || "_" || prediction_type AS score_name
        FROM peptide2prediction
    ''').fetchall()

    query = "CREATE VIEW v_predictions AS SELECT dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline, peptide.sequence, peptide.blosum62_2_similarity, variant.code, variant.variant_type, allele.allele_name"

    for i, (score_name,) in enumerate(scores_distinct):
        query += (
            f", MAX(CASE "
            f"WHEN peptide2prediction.prediction_tool || '_' || peptide2prediction.prediction_type = '{score_name}' "
            f"THEN peptide2prediction.value END) AS \"{score_name}\""
        )

    query += '''
        FROM peptide2prediction
        INNER JOIN peptide on peptide2prediction.peptide_id = peptide.id
        INNER JOIN allele on peptide2prediction.allele_id = allele.id
        INNER JOIN dataset on allele.dataset_id = dataset.id
        INNER JOIN variant on peptide.variant_id = variant.id
        GROUP BY dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline, peptide.sequence, peptide.blosum62_2_similarity, variant.code, variant.variant_type, allele.allele_name;
    '''
    cur.execute(query)

    ### View for gene signatures - v_gene_signatures
    cur.execute('DROP VIEW IF EXISTS v_gene_signatures')
    
    signatures_distinct = cur.execute(r'''
        SELECT DISTINCT gene_signature_type
        FROM gene_signature_value
    ''').fetchall()

    query = "CREATE VIEW v_gene_signatures AS SELECT dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline"

    for i, (sig_name,) in enumerate(signatures_distinct):
        query += (
            f", SUM(CASE "
            f"WHEN gene_signature_value.gene_signature_type = '{sig_name}' "
            f"THEN gene_signature_value.value ELSE 0 END) AS \"{sig_name}\""
        )

    query += '''
        FROM gene_signature_value
        INNER JOIN dataset on gene_signature_value.dataset_id = dataset.id
        GROUP BY dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline;
    '''

    cur.execute(query)
    
    ### View for variants
    variant_types = ["SNV", "SV", "fusion"]
    for variant in variant_types:
        cur.execute(f'DROP VIEW IF EXISTS v_variant_info_{variant.lower()}')
    
        variant_info_distinct = cur.execute(r'''
            SELECT DISTINCT variant_info.attribute_name
            FROM variant_info
            INNER JOIN variant on variant_info.variant_id = variant.id
            WHERE variant.variant_type = ?
        ''', (variant,)).fetchall()
        
        query = f"CREATE VIEW v_variant_info_{variant.lower()} AS SELECT dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline, variant.code, variant.variant_type"

        for i, (attr_name,) in enumerate(variant_info_distinct):
            query += (
                f", MAX(CASE "
                f"WHEN variant_info.attribute_name = '{attr_name}' "
                f"THEN variant_info.attribute_value ELSE NULL END) AS \"{attr_name}\""
            )
    
        query += f'''
            FROM variant_info
            INNER JOIN variant on variant_info.variant_id = variant.id
            INNER JOIN dataset on variant.dataset_id = dataset.id
            WHERE variant.variant_type = '{variant}'
            GROUP BY dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline, variant.code, variant.variant_type;
        '''
    
        cur.execute(query)
        
    ### View for MS prediction scores - v_predictions_ms
    cur.execute('DROP VIEW IF EXISTS v_predictions_ms')
    
    scores_distinct = cur.execute(r'''
        SELECT DISTINCT prediction_tool || "_" || prediction_type AS score_name
        FROM peptide_ms2prediction
    ''').fetchall()

    query = "CREATE VIEW v_predictions_ms AS SELECT dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline, peptide_ms.code, peptide_ms.peptide_type, peptide_ms2prediction.nmer_sequence, allele.allele_name"

    for i, (score_name,) in enumerate(scores_distinct):
        query += (
            f", MAX(CASE "
            f"WHEN peptide_ms2prediction.prediction_tool || '_' || peptide_ms2prediction.prediction_type = '{score_name}' "
            f"THEN peptide_ms2prediction.value END) AS \"{score_name}\""
        )

    query += '''
        FROM peptide_ms2prediction
        INNER JOIN peptide_ms on peptide_ms2prediction.peptide_id = peptide_ms.id
        INNER JOIN allele on peptide_ms2prediction.allele_id = allele.id
        INNER JOIN dataset on allele.dataset_id = dataset.id
        GROUP BY dataset.subject_id, dataset.glass_nr, dataset.treatment_timeline, peptide_ms.code, peptide_ms.peptide_type, peptide_ms2prediction.nmer_sequence, allele.allele_name;
    '''
    cur.execute(query)

def main():
    parser = argparse.ArgumentParser(description="Insert datasets and gene signatures into the database")
    parser.add_argument("--database_path", required=True, help="Path to the SQLite database file")

    args = parser.parse_args()

    con = sqlite3.connect(args.database_path)
    
    try:
        create_views(con)
        con.commit()
        
    except Exception as e:
        print(f"An error occurred when creating views: {e}. Rolling back all changes.")
        con.rollback()
        
    finally:
        con.close()
        
if __name__ == "__main__":
    main()
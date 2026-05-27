def calculate_similarity(gene, similarity = 0.8, min_similarity = 0.5, step = 0.1):
    '''
    This function calculates the functional similarity between a query gene and all 
    human genes in the database, based on shared GO terms annotated with the 
    'involved_in' relationship.
    
    The similarity is calculated as the proportion of the query gene's GO terms 
    that are shared with each candidate gene. If no results are found at the 
    requested similarity threshold, the threshold is automatically reduced by 
    the given step size until either results are found or the minimum threshold 
    is reached.

    Parameters:
        gene (str) -> The official gene symbol of the query gene.
        similarity (float) -> The initial similarity threshold to search for, 
                              expressed as a value between 0 and 1 (default = 0.8).
        min_similarity (float) -> The minimum similarity threshold accepted before 
                                  stopping the search (default: 0.5).
        step (float) -> The amount by which the similarity threshold is 
                        reduced at each iteration if no results are found 
                        (default: 0.1).

    Return:
        result (list of dict) -> A list of dictionaries, sorted by similarity in 
                                 descending order, where each dictionary contains:
                                    - gene_name (str): the candidate gene symbol.
                                    - common_terms (list): GO terms shared with the query gene.
                                    - n_common (int): number of shared GO terms.
                                    - similarity_pct (float): percentage of similarity 
                                      with respect to the query gene's GO terms.
                                    - all_terms (list): all GO terms of the candidate gene. 
    '''

    # Global items
    from pymongo import MongoClient
    RED='\033[91m'
    RESET='\033[0m'

    # Connect to Mongodb
    client = MongoClient('mongodb://localhost:27017/')
    db = client['saRNAs']
    collection = db['go_annotations']

    # Get the GO annotations for the target gene
    gene1 = collection.find_one({'gene_name':gene})

    if not gene1:
        print(f"Gene {gene} is not available in this database")
        client.close()
        return []
    
    go_terms = gene1['go_id']
    n_terms = len(go_terms)
    current_simi = similarity
    result = []

    while current_simi >= min_similarity:

        min_macth = round(n_terms * current_simi)
        print('*'*55)
        print(f"Target gene: {gene1['gene_name']}")
        print(f"Associated GO terms: {go_terms}")
        print(f"Minimun of coincidences: {min_macth}")
        print(f'Finding threshold: {current_simi*100}%')


        # Find genes that has unless 1 of this terms
        candidates = collection.find({
            'gene_name': {'$ne':gene},
            'go_id': {'$in': go_terms}
        })

        # Candidates that have the expected similarity
        result = []
        for candidate in candidates: 
            all_terms = candidate['go_id']
            c_terms = set(go_terms) & set(candidate['go_id'])
            n_common = len(c_terms)
            simi = n_common / n_terms

            if n_common >= min_macth:
                result.append({
                    'gene_name': candidate['gene_name'],
                    'common_terms': list(c_terms),
                    'n_common': n_common,
                    'similarity_pct': round(simi*100, 2),
                    'all_terms': all_terms
                })
        
        if result:
            if current_simi < similarity:
                print(f'\nNo results were found with similarity value of {similarity*100}%')
                print(f'There are results with a similarity value of {current_simi*100}%')
            break

        print(f'No results were found, modifying similarity...')
        print('*'*55)
        current_simi = round(current_simi - step, 1)

    if not result:
        print(f'No results were found even with the minimum similarity value of {min_similarity*100}')
        client.close()
        return[]
    
    # Sort similarity results
    result.sort(key=lambda x: x['similarity_pct'], reverse=True)

    client.close()
    print(f'\nNumber of genes found:{len(result)}')
    if n_terms <=2:
        print(f'{RED}WARNING: {gene} only has {n_terms} associated GO term(s)')
        print(f'Similarity results may not be biologically representative. Interpret them carefully.{RESET}')
    for g in result:
        print(g['gene_name'], g['similarity_pct'])
    return result



        

    

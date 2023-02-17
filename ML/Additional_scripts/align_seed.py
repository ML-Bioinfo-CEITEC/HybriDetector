import sys
from Bio.Seq import Seq
import re
import pandas as pd

def perfect_revcomp_match(driver, target):
    """
    Returns if the perfect reverse complement of the driver is in the target.
    """
    driver_revcomp = driver.reverse_complement()
    return target.find(driver_revcomp) != -1

def imperfect_revcomp_match(driver, target):
    """
    Returns if the reverse complmenet of the driver is imperfectly in the target and what positions had to be mutated to make a match.
    """
    mm_loc = [0]*len(driver)
    driver_revcomp = driver.reverse_complement()

    # we have an exact match
    if target.find(driver_revcomp) != -1:
        return (True, mm_loc)

    found = 0
    # try accepting mismatch separately at each position
    for i in range(len(driver)):
        mm_seq = re.compile(str(driver_revcomp[:i] + '.' + driver_revcomp[i+1:]))

        # TODO: findall can return multiple matches, but we check here for any match, not exact number
        if re.findall(mm_seq, str(target)):
            found += 1
            mm_loc[i] = 1

    if found > 0:
        # TODO: are we normalizing mismatch score so the sum is always 1?
        mm_loc = [loc/found for loc in mm_loc]
        found = True
    else:
        found = False

    return (found, mm_loc)

def count_prec_sens(df, metric):
    """
    Returns the precision and sensitivity of a given metric.
    """
    tp = len(df[((df[metric] == True) & (df['label'] == 1))])
    fp = len(df[((df[metric] == True) & (df['label'] == 0))])
    fn = len(df[((df[metric] == False) & (df['label'] == 1))])

    precision = (tp / (tp + fp))
    sensitivity = (tp / (tp + fn))

    return precision, sensitivity

# usage: python align_seed.py <path_to_input_clash_data> <path_to_output_metrics_csv_file> <path_to_precision_sensitivity_csv_file>
if __name__ == '__main__':

    if sys.argc != 4:
        print("Not enough arguments provided.")
        print("Usage: python align_seed.py <path_to_input_clash_data> <path_to_output_metrics_csv_file> <path_to_precision_sensitivity_csv_file>")

    # compute metrics for the data
    df = pd.read_csv(sys.argv[1], sep='\t')
    df['noncodingRNA'] = df['noncodingRNA'].apply(lambda x: Seq(x))
    df['gene'] = df['gene'].apply(lambda x: Seq(x))
    df['firstnt'] = df['noncodingRNA'].apply(lambda x: x[0])
    df['length'] = df['noncodingRNA'].apply(lambda x: len(x))
    df['seed'] = df.apply(lambda x: perfect_revcomp_match(x['noncodingRNA'][1:7], x['gene']), axis=1)
    df['stem9'] = df.apply(lambda x: perfect_revcomp_match(x['noncodingRNA'][0:9], x['gene']), axis=1)
    df['seed_1mm'] = df.apply(lambda x: (imperfect_revcomp_match(x['noncodingRNA'][1:7], x['gene']))[0], axis=1)
    df['stem9_1mm'] = df.apply(lambda x: (imperfect_revcomp_match(x['noncodingRNA'][0:9], x['gene']))[0], axis=1)
    df['maxL'] = 0
    for mirL in range(4, 15):
        df['seed1mm_' + str(mirL)] = df.apply(lambda x: (imperfect_revcomp_match(x['noncodingRNA'][1:mirL], x['gene']))[0], axis=1)
        df['maxL'] = df.apply(lambda x: mirL if (imperfect_revcomp_match(x['noncodingRNA'][1:mirL], x['gene']))[0] else x['maxL'], axis=1)
    df.to_csv(sys.argv[2], index=False)

    # compute the precision and sensitivity of the different metrics
    rows = []
    for metric in ['seed', 'stem9', 'seed_1mm', 'stem9_1mm'] + ['seed1mm_' + str(i) for i in range(4, 15)]:
        precision, sensitivity = count_prec_sens(df, metric)
        rows.append([metric, precision, sensitivity])
    df_metrics = pd.DataFrame(rows, columns=['metric', 'precision', 'sensitivity'])
    df_metrics.to_csv(sys.argv[3], index=False)

    # print out some statistics about the data
    print(df_metrics)

    print("FIRST NT:")
    print(
       pd.DataFrame({
            'firstnt_counts': df['firstnt'].value_counts().to_frame()['firstnt'], 
            'firstnt_percentage': df['firstnt'].value_counts(normalize=True).to_frame()['firstnt']})
    )

    print("LENGTH:")
    print(
       pd.DataFrame({
            'length_counts': df['length'].value_counts().to_frame()['length'], 
            'length_percentage': df['length'].value_counts(normalize=True).to_frame()['length']})
    )

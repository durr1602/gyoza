from snakemake.script import snakemake
import pandas as pd

def generate_read_stats(cutadapt_logfiles, pandaseq_logfiles, vsearch_logfiles, outpath):
    
    stats_dict = {}
    
    # Step 1 - Process cutadapt logs
    
    for f in cutadapt_logfiles:
        
        # Retrieve sample name
        sample_name = f[:-6].split('sample=')[1]
        
        # Parse cutadapt stats
        with open(f, 'r') as file:
            lines = file.readlines()
    
        # Check that log is properly formatted
        if 'This is cutadapt 5.0' not in lines[0]:
            raise Exception(f"Error.. {f} is not properly formatted. Make sure you've added the --use-conda flag in the snakemake command line, which specifies the correct package versions to be used")
        
        total_reads = int(lines[6].split('Total read pairs processed:')[1].strip().replace(',',''))
        r1_with_adapter = int(lines[7].split('Read 1 with adapter:')[1].split('(')[0].strip().replace(',',''))
        r2_with_adapter = int(lines[8].split('Read 2 with adapter:')[1].split('(')[0].strip().replace(',',''))
        trimmed_reads = int(lines[12].split('Pairs written (passing filters):')[1].split('(')[0].strip().replace(',',''))
        
        # Store all variables (number of reads) related to trimming step
        stats_dict[sample_name] = [total_reads, r1_with_adapter, r2_with_adapter, trimmed_reads]
    
    # Create dataframe to store all stats
    fullstats = pd.DataFrame.from_dict(stats_dict, orient='index',
                                       columns=['Total_raw_reads', 'R1_reads_with_adapter', 'R2_reads_with_adapter', 'Total_trimmed_reads'])
    
    # Calculate percentage of trimmed reads
    fullstats['Trimmed_%'] = fullstats['Total_trimmed_reads'] / fullstats['Total_raw_reads']
    
    # Step 2 - Process pandaseq logs
    
    for f in pandaseq_logfiles:
        
        # Retrieve sample name
        sample_name = f[:-6].split('sample=')[1]
        
        # Check that log is properly formatted
        with open(f, 'r') as file:
            first_line = file.readline()
        if 'INFO	VER	pandaseq 2.11' not in first_line:
            raise Exception(f"Error.. {f} is not properly formatted. Make sure you've added the --use-conda flag in the snakemake command line, which specifies the correct package versions to be used")
        
        # Parse pandaseq stats
        logfile = pd.read_csv(f, sep='\t', skiprows=25, skipfooter=1, engine='python', names=['id','err_stat','field','value','details'])
        stats = logfile[logfile.field.isin(['LOWQ','NOALGN','OK','READS','SLOW'])].iloc[-5:,:][['field','value']]
        stats['value'] = stats.value.astype(int)
        
        # Validation step - check that the number of processed reads corresponds from trim output to merge input
        if fullstats.loc[fullstats.index == sample_name, 'Total_trimmed_reads'].item() != stats.loc[stats.field == 'READS', 'value'].item():
            print('---ERROR---\nNumber of written reads in trim output does not correspond to number of processed reads in merge input')
        
        # Add stats to dataframe
        fullstats.loc[fullstats.index == sample_name,
                      ['Not_merged_LOWQ','Not_merged_NOALGN','Total_merged_reads']] = stats.set_index('field').T[['LOWQ','NOALGN','OK']].values
    
    # Cast type to remove ','
    fullstats = fullstats.astype({'Not_merged_LOWQ': 'int',
                                  'Not_merged_NOALGN': 'int',
                                  'Total_merged_reads': 'int'
                                 })
    
    # Calculate percentage of properly merged reads relative to number of trimmed reads
    fullstats['Merged_%'] = fullstats['Total_merged_reads'] / fullstats['Total_trimmed_reads']
    
    # Step 3 - Process vsearch logs
    
    for f in vsearch_logfiles:
    
        # Retrieve sample name
        sample_name = f[:-6].split('sample=')[1]
        
        # Parse vsearch stats
        with open(f, 'r') as file:
            lines = file.readlines()

        # Check that log is properly formatted
        if 'vsearch v2.29.3' not in lines[0]:
            raise Exception(f"Error.. {f} is not properly formatted. Make sure you've added the --use-conda flag in the snakemake command line, which specifies the correct package versions to be used")        
        
        # Add singletons total   
        singletons = lines[-1].split(' clusters discarded')[0].split(',')[-1].strip()
        fullstats.loc[fullstats.index == sample_name, 'Nb_singletons'] = singletons
    
    # Convert column type to int
    fullstats['Nb_singletons'] = fullstats['Nb_singletons'].astype(int)
    
    # Calculate number of reads corresponding to aggregated sequences    
    fullstats['Nb_non-singletons'] = fullstats['Total_merged_reads'] - fullstats['Nb_singletons']
    
    # Calculate percentage of reads corresponding to non-singletons (relative to number of properly merged reads)
    fullstats['Aggregated_%'] = fullstats['Nb_non-singletons'] / fullstats['Total_merged_reads']
    
    # Output to csv file
    fullstats.to_csv(outpath)
    
    return

generate_read_stats(snakemake.input.cutadapt_logs,
                    snakemake.input.pandaseq_logs,
                    snakemake.input.vsearch_logs,
                    snakemake.output[0])

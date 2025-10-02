// dedicated to the utility functions
class ViroSeekUtils {

    // Function to chek if we
    public static Boolean is_url(file) {
        
        // check if the file is a URL
        if (file =~ /^(http|https|ftp|s3|az|gs):\/\/.*/) {
            return true
        } else {
            return false
        }
    }
    // Function to chek if we
    public static Boolean is_fastq(file) {
        
        // check if the file is a URL
        if (file =~ /.*(fq|fastq|fq.gz|fastq.gz)$/) {
            return true
        } else {
            return false
        }
    }

    // Function to extract the basename of a file
    public static String getCleanName(file) {
        def fileClean = file[0].baseName.replaceAll(/\.(gz)$/, '') // remove .gz
        fileClean = fileClean.replaceAll(/\.(fasta|fa)$/, '') // remove .fasta or .fa
        fileClean = fileClean.replaceAll(/\.(fastq|fq)$/, '') // remove .fastq or .fq
        return fileClean
    }

    /*
    * Validate a parameter (local file path or URL).
    * Returns a Nextflow `file()` handle, or null if the parameter is not provided.
    *
    * @param paramValue  the parameter value (URL or local path)
    * @param label       human-readable name used in error messages
    * @return            a Nextflow file() handle, or null when paramValue is null/empty
    */
    static def resolveFile(paramValue, String label) {
        if (!paramValue) return null

        if (is_url(paramValue)) {
            ref_file = file(paramValue)
            ref_file = ref_file.getName()
        } else {
            def f = new File(paramValue)
            if (!f.exists()) {
                exit 1, "Error: ${label} file <${paramValue}> does not exist.\n"
            }
            ref_file = f.getName()
            return ref_file
        }
    }
}
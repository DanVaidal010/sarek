/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOG INICIAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// Para obtener información para la fecha y el timestamp:
import java.text.SimpleDateFormat
def date = new Date()
def sdf  = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss")

// Parámetros para el log:
params.date = sdf.format(date)

// Información de log mostrada en pantalla:
log.info """
    S T A T I S T I C S  W O R K F L O W
    =====================================
    Nextflow version  = ${nextflow.version}
    Start timestamp   = ${params.date}
    Project directory = ${projectDir}
    Output directory  = ${params.outdir}
    """
    .stripIndent()

// Información de ayuda:
params.help = false

if (params.help) { // Se muestra en pantalla cuando se introduce la flag "--help".
    log.info "Texto de ayuda."
    exit 1
}

workflow.onComplete {
    def end_date = new Date()
    params.end_date = sdf.format(end_date)
    log.info "\nCompletion timestamp: ${params.end_date} " // Better than ${workflow.complete}
    log.info(workflow.success ? "\nLa ejecución del pipeline ha terminado correctamente.\n" : "Oops, algo ha fallado.")
}

workflow.onError {
    log.info "\nOops... La ejecución del pipeline se ha detenido con el siguiente mensaje: ${workflow.errorMessage}"
}

include { STATISTICS_SUBWORKFLOW } from '/media/scratch0/20230622_UB110_pruebas_sarek/02estadisticas/nextflow_statistics_workflow/statistics_subworkflow.nf'

workflow {
    intervals_capture = params.intervals_capture  ? Channel.value(params.intervals_capture)      : null
    intervals_target  = params.intervals_target   ? Channel.value(params.intervals_target)       : null
    sarek_directory   = params.sarek_analysis_dir ? Channel.value(params.sarek_analysis_dir)     : null
    fasta             = params.fasta              ? Channel.fromPath(params.fasta).collect()     : null
    fai               = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect() : null

    csv_file = params.samplesheet
    samples  = extract_csv(file(csv_file))

    STATISTICS_SUBWORKFLOW(samples, intervals_capture, intervals_target, sarek_directory, fasta, fai)
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
def extract_csv(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, samplesheet_line_count = 0;
        while ((line = reader.readLine()) != null) {samplesheet_line_count ++}
            if (samplesheet_line_count < 2) {
                error("Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines.")
            }
    }

    // Additional check of sample sheet:
    // 1. Each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def patient_sample_lane_combinations = []
    def sample2patient = [:]

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
            if ( !row.lane ) {  // This also handles the case where the lane is left as an empty string
                error('The sample sheet should specify a lane for patient "' + row.patient.toString() + '" and sample "' + row.sample.toString() + '".')
            }
            def patient_sample_lane = [row.patient.toString(), row.sample.toString(), row.lane.toString()]
            if (patient_sample_lane in patient_sample_lane_combinations) {
                error('The patient-sample-lane combination "' + row.patient.toString() + '", "' + row.sample.toString() + '", and "' + row.lane.toString() + '" is present multiple times in the sample sheet.')
            } else {
                patient_sample_lane_combinations.add(patient_sample_lane)
            }

            if (!sample2patient.containsKey(row.sample.toString())) {
                sample2patient[row.sample.toString()] = row.patient.toString()
            } else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
                error('The sample "' + row.sample.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sample.toString()] + '" in the sample sheet.')
            }
        }

    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0

    Channel.of(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all ++
            if (!(row.patient && row.sample)) {
                error("Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'.")
            }
            else if (row.patient.contains(" ") || row.sample.contains(" ")) {
                error("Invalid value in csv file. Values for 'patient' and 'sample' can not contain space.")
            }
            [ [ row.patient.toString(), row.sample.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map{ row, num_lanes -> // from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 1) sample_count_tumor ++
        else sample_count_normal ++

        // start from BAM
        if (row.lane && row.bam) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def bai         = row.bai ? file(row.bai,   checkIfExists: true) : []

            meta.num_lanes  = num_lanes.toInteger()
            meta.data_type  = 'bam'

            meta.size       = 1 // default number of splitted fastq

            return [ meta, bam, bai ]
        }

        // start from fastq
        if (row.lane && row.fastq_2) {
            // meta.id         = "${row.sample}-${row.lane}".toString()
            meta.id         = row.sample.toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            // def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''

            // def flowcell    = flowcellLaneFromFastq(fastq_1)
            // Don't use a random element for ID, it breaks resuming
            // def read_group  = "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            meta.num_lanes  = num_lanes.toInteger()
            // meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'

            meta.size       = 1 // default number of splitted fastq

            return [ meta, fastq_1, fastq_2 ]
        }
    }
}

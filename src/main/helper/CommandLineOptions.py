class CommandLineOptions:
    path_to_fasta_sequence_file = None
    path_to_sequence_file2 = None
    path_to_substitution_matrix = None
    open_gap_cost = None
    extend_gap_cost = None
    traceback_complete = None
    path_to_output = None

    def setAlignmentOptions(self,
                                    fasta_file1,
                                    fasta_file2,
                                    substitution_matrix,
                                    open_gap_cost,
                                    extend_gap_cost,
                                    traceback_complete,
                                    outputfile):
        self.path_to_fasta_sequence_file = fasta_file1
        self.path_to_sequence_file2 = fasta_file2
        self.extend_gap_cost = extend_gap_cost
        self.open_gap_cost = open_gap_cost
        self.path_to_substitution_matrix = substitution_matrix
        self.path_to_output = outputfile
        self.traceback_complete = traceback_complete

    def setRNAFoldingOptions(self, inputfile, outputfile):
        self.path_to_sequence_file2 = inputfile
        self.path_to_output = outputfile


class MultiAlignmentResult:
    id = None
    score = None
    aligned_sequence_list = None

    def __init__(self):
        self.id = []
        self.aligned_sequence_list = []

    def alignment_result (self, id, align_sequence_list, score):
        self.id = id
        self.aligned_sequence_list = align_sequence_list
        self.score =score

    def to_tuple(self, **sequences):
        for i in range (len(**sequences)):
            self.id.append((i, sequences(i)))
        return (self.id, self.score, self.aligned_sequence_list)
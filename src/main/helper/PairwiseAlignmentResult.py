class PairwiseAlignmentResult:
    id1 = None
    id2 = None
    score = None
    aligned_sequence_list = None
    aligned_sequence1 =None
    aligned_sequence2= None
    def __init__(self):
        self.aligned_sequence1 =[]
        self.aligned_sequence2 =[]
        self.aligned_sequence_list = [self.aligned_sequence1, self.aligned_sequence2]

    def alignment_result_list(self, id1, id2, score, aligned_sequence_list):
        self.id1= id1
        self.id2 = id2
        self.score = score
        self.aligned_sequence_list = aligned_sequence_list
        if len(self.aligned_sequence_list) > 0:
            self.aligned_sequence1 = self.aligned_sequence_list[0][0]
            self.aligned_sequence2 = self.aligned_sequence_list[0][1]
        else:
            self.aligned_sequence1 = None
            self.aligned_sequence2 = None

    def alignment_result(self, id1, id2, score, aligned_sequence1, aligned_sequence2):
        self.id1 = id1
        self.id2 = id2
        self.score = score
        self.aligned_sequence1 = aligned_sequence1
        self.aligned_sequence2 = aligned_sequence2
        self.aligned_sequence_list = [[aligned_sequence1, aligned_sequence2]]

    def to_tuple(self, seq1, seq2):
        return (self.id1, seq1, self.id2, seq2, self.score, self.aligned_sequence_list)


class PairwiseAlignmentResult:
    id1 = None
    id2 = None
    score = None
    aligned_sequence_list = None
    aligned_sequence1 = None
    aligned_sequence2 = None

    def __init__(self):
        self.aligned_sequence1 = []
        self.aligned_sequence2 = []
        self.aligned_sequence_list = [self.aligned_sequence1, self.aligned_sequence2]

    def alignment_result_list(self, id1, id2, score, aligned_sequence_list):
        self.id1 = id1
        self.id2 = id2
        self.score = score
        self.aligned_sequence_list = aligned_sequence_list
        if len(self.aligned_sequence_list) > 0:
            self.aligned_sequence1 = self.aligned_sequence_list[0][0]
            self.aligned_sequence2 = self.aligned_sequence_list[0][1]
        else:
            self.aligned_sequence1 = None
            self.aligned_sequence2 = None

    def alignment_result(self, id1, id2, score, aligned_sequence1, aligned_sequence2):
        self.id1 = id1
        self.id2 = id2
        self.score = score
        self.aligned_sequence1 = aligned_sequence1
        self.aligned_sequence2 = aligned_sequence2
        self.aligned_sequence_list = [[aligned_sequence1, aligned_sequence2]]

    def to_tuple(self, seq1, seq2):
        return (self.id1, seq1, self.id2, seq2, self.score, self.aligned_sequence_list)

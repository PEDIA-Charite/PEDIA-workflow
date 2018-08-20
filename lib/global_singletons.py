from lib.api import mutalyzer, jannovar, omim, phenomizer, aws_download
from lib import errorfixer

MUTALYZER_INST = mutalyzer.Mutalyzer()

JANNOVAR_INST = jannovar.JannovarClient()

OMIM_INST = omim.Omim()

PHENOMIZER_INST = phenomizer.PhenomizerService()

ERRORFIXER_INST = errorfixer.ErrorFixer()

AWS_INST = aws_download.AWSBucket()

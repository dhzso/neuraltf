#!/usr/bin/env python
import sys
from pptx import Presentation
ppt_path = r'D:/Bioinformatics/projects/NeuralTF/docs/NeuralTF_SOTA_Final.pptx'
prs = Presentation(ppt_path)
print('Slide count:', len(prs.slides))

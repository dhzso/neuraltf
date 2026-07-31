"""Build v4-v6 bridge CSV for TF genes.

Usage: python scripts/build_bridge.py
"""
exec(compile(open(r"D:\Bioinformatics\src\bioforge\projects\neuraltf\pipeline.py").read() +
    "\nNeuralTFPipeline().load_reference_tables()", "build_bridge", "exec"))

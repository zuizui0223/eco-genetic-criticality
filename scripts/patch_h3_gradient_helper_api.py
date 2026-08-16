from pathlib import Path

chain = Path("causal_model/mutation_primary_h1_h2_h3_chain.py")
text = chain.read_text(encoding="utf-8")
old = """    canonical = canonical_h1_certificate(calibration.parameters.interaction_feedback, one_large.patch_areas[0], calibration.parameters.area_reference,\n        (interval[0] + interval[1]) / 2.0, base)\n"""
new = """    canonical = canonical_h1_certificate(\n        feedback_strength=calibration.parameters.interaction_feedback,\n        area=one_large.patch_areas[0],\n        area_reference=calibration.parameters.area_reference,\n        barrier=(interval[0] + interval[1]) / 2.0,\n        trait_parameters=base,\n    )\n"""
if old not in text:
    if new in text:
        print("helper already uses keyword-only canonical API")
    else:
        raise SystemExit("expected canonical_h1_certificate call not found")
else:
    chain.write_text(text.replace(old, new, 1), encoding="utf-8")
    print("patched mutation-primary helper to keyword-only canonical API")

# Audio safety math: Python bindings

Python bindings for [audio-safety-math](https://github.com/chethandvg/audio-safety-math), built from the Rust core with PyO3 and maturin. Function names and results are identical to the Rust crate; the shared golden vectors check it in CI. Signals and lists of levels are accepted as numpy arrays or any sequence of numbers, and signals come back as numpy arrays. The dose functions default to the NIOSH rule (85 dBA, 3 dB exchange rate, 8 hours).

```bash
pip install "git+https://github.com/chethandvg/audio-safety-math@v0.1.0#subdirectory=bindings/python"
```

```python
import audio_safety_math as asm

asm.allowed_time(88.0)  # 4.0 hours
asm.sum_incoherent_db([85.0, 85.0])  # 88.0103 dB
```

These are formulas and reference values, not medical advice.

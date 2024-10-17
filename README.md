# Constrained Dynamics Simulation

## Sources
https://www.cs.cmu.edu/~baraff/pbm/constraints.pdf
https://en.wikipedia.org/wiki/Minimal_residual_method
https://en.wikipedia.org/wiki/Conjugate_gradient_squared_method

## Installation/Usage (Linux & macOS)

```sh
git clone https://github.com/philipobi/constrained-dynamics.git
cd constrained-dynamics
mkdir bin
make build
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib

./bin/main > positions.txt
python3 animation.py
```

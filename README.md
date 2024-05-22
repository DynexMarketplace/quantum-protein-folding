# Quantum Protein Folding
Protein folding is the physical process by which a protein chain attains its functional three-dimensional structure from a simple sequence of amino acids. This folding occurs spontaneously, guided by interactions among the amino acids and the surrounding environment, which determine the protein's final shape. Correct folding is crucial for a protein's function, as misfolding can lead to diseases like Alzheimer's and Parkinson's.

# Usage

1. Launch in Github Codespaces and wait until the codepsace is fully initialised

2. Add your account keys by drag&drop of your dynex.ini into the main folder

3. Verify your account keys by typing the following command in the console:

```
python
>>> import dynex
>>> dynex.test()
>>> exit()
```

Your console will perform tests and validate your account keys. You should see the following message:

```
[DYNEX] TEST RESULT: ALL TESTS PASSED
```

4. Run the demo by typing the following command:

```
python main.py
```

The program will output and save the folded protein in the file "result.png".



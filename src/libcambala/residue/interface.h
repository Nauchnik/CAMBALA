#ifndef RESIDUE_CALCULATOR_H_
#define RESIDUE_CALCULATOR_H_

class ResidueCalculator
{
public:
	virtual void LoadImmutableData () = 0;
	virtual void CalculatePointResidue() = 0;
	virtual ~ResidueCalculator(){};
};

#endif

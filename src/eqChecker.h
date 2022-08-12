#ifndef _PEQCHECKER_H_
#define _PEQCHECKER_H_

#include "bddSystem.h"

enum class EqType 
{ 
    Feq, // Full equivalence 
    Peq, // Partial equivalence 
    PeqS // Partial equivalence Special
};

class EquivalenceChecker : public BDDSystem
{
public:
    // Constructor and Destructor
    EquivalenceChecker
    (
        std::vector<std::vector<GateType>>& gates,
        std::vector<std::vector<std::vector<int>>>& qubits,
        int n,
        int nQin, 
        int nQout, 
        int r, 
        bool isReorder,
        EqType eqType
    );

    ~EquivalenceChecker()  
    {
        clear();
    }

    void check();
    void printInfo(double runtime, size_t memPeak) const;

private:
    std::vector<std::vector<GateType> > _gates;             // gates in circuits. [nCircuit]*[#gate]
    std::vector<std::vector<std::vector<int> > > _qubits;   // ith qubits of gates in circuits. [nCircuit]*[#gates]*[#qubits]
    int _nQin;                                              // #data qubits (d in the paper)
    int _nQout;                                             // #measured qubits (m in the paper)     
    int _ratio;                                             // gate count ratio. |circuit2|/|circuit1|
    bool _isEq;                                             // if the result is equivalent or not.
    EqType _eqType;                                         // the equivalence type (Feq/Peq/PeqS)  

    void invertCircuit(std::vector<GateType> &gate);
    void init();
    void applyGate(int ithCircuit, GateType type, std::vector<int> qubit, bool right);
    void calculateMiter();
    void extract(int ithCircuit);
    void checkFeq();
    void checkPeq();
    void checkPeqS();
    void printResult() const;

    // Clean up EquivalenceChecker
    void clear() 
    {
        _gates[0].clear();
        _gates[1].clear();
        _gates.clear();
        _qubits[0].clear();
        _qubits[1].clear();
        _qubits.clear();
    };
};

#endif

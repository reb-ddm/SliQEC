OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
cx q[5], q[9];
cx q[8], q[5];
t q[9];
s q[1];
t q[4];
s q[0];
s q[5];
t q[7];
s q[2];
h q[8];
ccx q[8], q[1], q[2];
cx q[6], q[9];
h q[7];
ccx q[4], q[6], q[1];
s q[3];
ccx q[5], q[6], q[2];
ccx q[4], q[5], q[6];
ccx q[6], q[4], q[5];
ccx q[0], q[6], q[9];
t q[9];
cx q[6], q[3];
s q[1];
s q[6];
s q[4];
ccx q[8], q[7], q[3];
t q[4];
s q[3];
ccx q[5], q[3], q[0];
ccx q[4], q[9], q[7];
ccx q[4], q[8], q[6];
tdg q[0];
y q[2];
sdg q[3];
s q[7];
cx q[9], q[8];
h q[8];
h q[7];
h q[8];

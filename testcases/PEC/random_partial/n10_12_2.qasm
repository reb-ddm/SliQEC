OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
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
h q[4];
t q[7];
s q[8];
h q[7];
ccx q[3], q[6], q[8];
cx q[2], q[8];
ccx q[6], q[1], q[2];
s q[7];
cx q[6], q[0];
t q[2];
cx q[8], q[7];
s q[1];
ccx q[8], q[4], q[7];
ccx q[6], q[4], q[2];
t q[8];
cx q[6], q[0];
t q[8];
h q[9];
ccx q[4], q[3], q[0];
cx q[8], q[1];
ccx q[9], q[5], q[0];
ccx q[9], q[2], q[3];
s q[6];
t q[4];
h q[3];
ccx q[8], q[4], q[1];
s q[9];
t q[2];
cx q[5], q[1];
t q[4];
cx q[0], q[1];
sdg q[0];
sdg q[1];
cx q[0], q[1];
sdg q[0];
sdg q[2];
sdg q[3];
cx q[2], q[3];
tdg q[2];
cx q[8], q[5];
h q[6];
cx q[6], q[5];
h q[6];
cx q[8], q[5];
cx q[10], q[7];
cx q[11], q[9];

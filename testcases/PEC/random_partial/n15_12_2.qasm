OPENQASM 2.0;
include "qelib1.inc";
qreg q[18];
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
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
ccx q[1], q[4], q[11];
s q[1];
cx q[0], q[14];
t q[14];
s q[14];
h q[6];
t q[7];
h q[8];
t q[12];
s q[5];
h q[1];
h q[7];
s q[9];
t q[6];
s q[0];
s q[3];
cx q[11], q[7];
t q[8];
ccx q[2], q[10], q[5];
s q[7];
h q[10];
ccx q[5], q[11], q[12];
ccx q[11], q[13], q[9];
h q[13];
h q[12];
t q[14];
h q[11];
cx q[10], q[8];
ccx q[14], q[11], q[6];
ccx q[8], q[4], q[12];
h q[1];
ccx q[1], q[9], q[6];
h q[8];
t q[13];
h q[0];
s q[0];
t q[0];
t q[7];
ccx q[10], q[13], q[6];
t q[9];
h q[12];
cx q[11], q[4];
cx q[8], q[7];
cx q[6], q[14];
ccx q[14], q[8], q[3];
tdg q[1];
sdg q[0];
cx q[0], q[1];
sdg q[2];
cx q[3], q[2];
sdg q[2];
cx q[3], q[2];
sdg q[2];
sdg q[4];
t q[9];
t q[11];
cx q[14], q[10];
cx q[7], q[12];
ccx q[8], q[13], q[12];
cx q[7], q[13];
ccx q[12], q[9], q[13];
cx q[7], q[14];
cx q[15], q[12];
cx q[16], q[6];
cx q[17], q[11];

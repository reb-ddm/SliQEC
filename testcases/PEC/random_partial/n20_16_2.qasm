OPENQASM 2.0;
include "qelib1.inc";
qreg q[23];
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
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
cx q[7], q[0];
ccx q[11], q[5], q[16];
t q[6];
ccx q[0], q[6], q[19];
h q[17];
h q[13];
h q[14];
h q[10];
cx q[18], q[13];
ccx q[14], q[10], q[19];
cx q[18], q[5];
ccx q[6], q[16], q[7];
s q[3];
ccx q[3], q[16], q[17];
h q[15];
ccx q[0], q[13], q[8];
cx q[17], q[7];
ccx q[9], q[15], q[6];
s q[13];
t q[11];
h q[17];
cx q[7], q[8];
cx q[5], q[13];
t q[5];
s q[18];
t q[16];
s q[14];
cx q[0], q[2];
h q[2];
s q[15];
t q[0];
ccx q[6], q[11], q[3];
h q[14];
s q[1];
s q[3];
s q[18];
h q[4];
h q[9];
ccx q[10], q[3], q[4];
t q[18];
s q[8];
t q[18];
t q[6];
h q[4];
t q[3];
s q[16];
ccx q[11], q[9], q[10];
cx q[14], q[6];
ccx q[3], q[0], q[6];
ccx q[5], q[18], q[10];
s q[17];
h q[4];
t q[6];
cx q[19], q[18];
cx q[15], q[6];
cx q[10], q[3];
h q[17];
h q[12];
s q[14];
h q[3];
cx q[2], q[1];
tdg q[2];
tdg q[1];
cx q[2], q[1];
sdg q[1];
sdg q[3];
tdg q[3];
cx q[5], q[6];
tdg q[5];
z q[7];
cx q[8], q[9];
tdg q[8];
t q[19];
t q[18];
s q[17];
t q[16];
h q[17];
h q[15];
cx q[18], q[12];
cx q[16], q[19];
h q[11];
t q[11];
cx q[20], q[16];
cx q[21], q[13];
cx q[22], q[11];

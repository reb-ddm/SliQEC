OPENQASM 2.0;
include "qelib1.inc";
qreg q[25];
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
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[15];
ccx q[15], q[21], q[0];
s q[12];
t q[0];
cx q[23], q[0];
h q[2];
ccx q[11], q[5], q[20];
s q[18];
cx q[1], q[2];
ccx q[0], q[21], q[20];
h q[5];
h q[21];
h q[14];
cx q[8], q[1];
ccx q[18], q[24], q[4];
h q[24];
s q[16];
h q[17];
cx q[8], q[12];
t q[1];
h q[2];
s q[20];
ccx q[13], q[16], q[5];
s q[24];
ccx q[20], q[9], q[6];
t q[19];
ccx q[24], q[11], q[23];
ccx q[0], q[12], q[14];
h q[6];
s q[14];
ccx q[11], q[22], q[6];
ccx q[18], q[9], q[16];
ccx q[21], q[23], q[13];
t q[15];
ccx q[16], q[11], q[3];
h q[17];
s q[5];
t q[9];
cx q[18], q[14];
t q[0];
cx q[15], q[2];
s q[21];
ccx q[15], q[4], q[12];
ccx q[8], q[24], q[14];
h q[23];
s q[5];
cx q[14], q[24];
h q[1];
h q[11];
t q[10];
h q[6];
ccx q[14], q[7], q[19];
ccx q[10], q[8], q[14];
ccx q[12], q[21], q[1];
cx q[10], q[18];
h q[1];
t q[14];
h q[12];
cx q[7], q[1];
ccx q[7], q[13], q[1];
t q[6];
t q[13];
cx q[8], q[7];
s q[15];
ccx q[7], q[22], q[16];
h q[21];
t q[0];
s q[6];
s q[15];
s q[4];
h q[22];
cx q[8], q[11];
s q[22];
cx q[18], q[6];
s q[21];
z q[0];
tdg q[1];
cx q[2], q[1];
tdg q[1];
sdg q[1];
sdg q[3];
cx q[6], q[5];
sdg q[5];
tdg q[5];
cx q[6], q[5];
sdg q[5];
cx q[8], q[7];
sdg q[7];
cx q[8], q[7];
sdg q[7];
sdg q[9];
ccx q[23], q[15], q[16];
ccx q[13], q[17], q[20];
ccx q[13], q[22], q[18];
t q[21];
s q[21];
ccx q[24], q[20], q[16];
s q[23];
t q[19];
t q[23];
t q[21];
h q[14];
h q[22];
cx q[21], q[13];
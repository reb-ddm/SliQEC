OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
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
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
t q[29];
t q[10];
cx q[23], q[5];
cx q[21], q[4];
cx q[4], q[13];
ccx q[20], q[11], q[10];
h q[0];
s q[22];
t q[10];
s q[3];
t q[11];
ccx q[1], q[3], q[22];
cx q[27], q[20];
ccx q[27], q[20], q[19];
s q[0];
cx q[15], q[2];
s q[2];
s q[6];
ccx q[6], q[11], q[0];
t q[22];
t q[27];
h q[10];
ccx q[1], q[11], q[5];
cx q[2], q[23];
cx q[29], q[17];
s q[21];
cx q[15], q[17];
t q[0];
t q[0];
h q[14];
cx q[22], q[0];
s q[16];
h q[21];
h q[29];
cx q[3], q[13];
s q[27];
cx q[18], q[29];
t q[18];
cx q[0], q[27];
ccx q[1], q[18], q[6];
s q[19];
s q[14];
s q[14];
s q[0];
ccx q[24], q[1], q[29];
ccx q[15], q[28], q[12];
t q[26];
cx q[10], q[19];
ccx q[17], q[20], q[3];
s q[18];
h q[14];
cx q[28], q[15];
t q[13];
s q[24];
h q[1];
t q[12];
cx q[7], q[19];
ccx q[0], q[17], q[20];
s q[18];
t q[28];
s q[22];
s q[25];
t q[17];
ccx q[22], q[5], q[13];
ccx q[29], q[3], q[22];
t q[12];
s q[2];
ccx q[29], q[11], q[26];
cx q[4], q[24];
cx q[25], q[17];
cx q[28], q[3];
t q[13];
t q[25];
t q[28];
cx q[11], q[26];
cx q[21], q[8];
h q[6];
cx q[11], q[6];
h q[10];
h q[12];
cx q[6], q[18];
cx q[7], q[29];
cx q[10], q[3];
h q[6];
t q[2];
s q[9];
h q[6];
s q[7];
ccx q[22], q[24], q[13];
t q[9];

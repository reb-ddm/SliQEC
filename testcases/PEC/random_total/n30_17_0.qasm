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
s q[29];
ccx q[1], q[14], q[8];
t q[28];
ccx q[3], q[20], q[12];
t q[10];
cx q[10], q[12];
t q[26];
t q[6];
ccx q[29], q[9], q[21];
cx q[27], q[23];
s q[1];
h q[8];
t q[27];
s q[19];
ccx q[15], q[28], q[9];
cx q[7], q[27];
h q[25];
s q[29];
ccx q[8], q[5], q[29];
cx q[28], q[9];
ccx q[21], q[8], q[0];
t q[3];
ccx q[11], q[21], q[10];
ccx q[22], q[27], q[10];
ccx q[7], q[23], q[27];
ccx q[7], q[8], q[0];
ccx q[5], q[10], q[22];
ccx q[23], q[6], q[19];
t q[11];
s q[27];
cx q[11], q[27];
ccx q[12], q[27], q[10];
t q[21];
cx q[12], q[9];
s q[20];
ccx q[22], q[3], q[8];
h q[10];
cx q[29], q[27];
s q[28];
ccx q[15], q[17], q[27];
h q[8];
s q[5];
h q[8];
t q[17];
h q[10];
ccx q[0], q[4], q[11];
cx q[3], q[19];
s q[12];
s q[8];
h q[0];
ccx q[26], q[11], q[17];
ccx q[28], q[21], q[8];
t q[29];
h q[16];
ccx q[0], q[18], q[12];
s q[14];
cx q[9], q[28];
s q[12];
h q[6];
s q[13];
h q[25];
h q[21];
t q[2];
t q[25];
s q[13];
h q[3];
s q[3];
s q[25];
s q[6];
ccx q[11], q[27], q[2];
ccx q[7], q[21], q[25];
t q[24];
s q[16];
s q[8];
h q[6];
ccx q[22], q[25], q[27];
cx q[7], q[2];
ccx q[11], q[26], q[24];
h q[18];
t q[21];
s q[28];
s q[0];
ccx q[20], q[22], q[25];
cx q[24], q[23];
ccx q[20], q[16], q[11];
ccx q[26], q[25], q[7];
cx q[23], q[4];
ccx q[11], q[24], q[25];
ccx q[16], q[5], q[15];
h q[6];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[45];
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
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
ccx q[33], q[19], q[11];
cx q[32], q[22];
ccx q[5], q[28], q[34];
cx q[33], q[16];
t q[37];
h q[0];
h q[16];
t q[16];
cx q[17], q[34];
h q[2];
s q[13];
ccx q[28], q[12], q[8];
h q[2];
cx q[8], q[7];
s q[8];
h q[5];
cx q[26], q[35];
s q[6];
ccx q[21], q[28], q[12];
t q[19];
ccx q[20], q[30], q[8];
h q[22];
s q[20];
ccx q[23], q[11], q[31];
ccx q[38], q[35], q[32];
h q[30];
s q[24];
ccx q[14], q[33], q[3];
t q[4];
s q[25];
h q[7];
ccx q[10], q[31], q[12];
ccx q[13], q[28], q[39];
s q[38];
s q[2];
ccx q[7], q[37], q[36];
cx q[28], q[0];
ccx q[11], q[31], q[12];
t q[17];
h q[23];
cx q[0], q[4];
ccx q[22], q[16], q[34];
s q[37];
ccx q[8], q[0], q[38];
ccx q[27], q[36], q[14];
h q[31];
h q[33];
s q[17];
cx q[22], q[17];
s q[19];
cx q[4], q[8];
cx q[38], q[15];
s q[36];
s q[20];
ccx q[6], q[14], q[11];
cx q[25], q[31];
cx q[14], q[25];
s q[17];
h q[7];
ccx q[2], q[36], q[10];
cx q[11], q[19];
cx q[38], q[23];
s q[29];
ccx q[16], q[38], q[21];
cx q[29], q[33];
s q[35];
h q[21];
t q[0];
t q[30];
cx q[23], q[19];
t q[27];
t q[17];
h q[39];
t q[36];
t q[23];
t q[13];
ccx q[23], q[12], q[0];
cx q[6], q[24];
s q[39];
cx q[6], q[38];
s q[33];
h q[33];
t q[19];
cx q[2], q[17];
s q[32];
cx q[9], q[31];
t q[14];
s q[14];
ccx q[22], q[25], q[18];
ccx q[16], q[24], q[9];
s q[15];
t q[34];
cx q[29], q[20];
s q[8];
t q[4];
t q[4];
h q[2];
cx q[15], q[30];
h q[32];
h q[15];
ccx q[4], q[7], q[1];
h q[1];
t q[28];
s q[19];
h q[19];
cx q[26], q[38];
t q[35];
t q[20];
ccx q[10], q[36], q[25];
t q[39];
s q[7];
ccx q[23], q[16], q[37];
h q[29];
s q[35];
t q[8];
h q[25];
t q[14];
t q[33];
ccx q[14], q[15], q[7];
s q[1];
cx q[1], q[0];
sdg q[0];
cx q[1], q[0];
tdg q[0];
sdg q[2];
y q[4];
y q[5];
x q[6];
sdg q[7];
sdg q[10];
cx q[9], q[10];
sdg q[9];
sdg q[11];
cx q[12], q[11];
tdg q[11];
sdg q[14];
sdg q[13];
cx q[13], q[14];
sdg q[13];
cx q[16], q[15];
tdg q[15];
x q[15];
cx q[16], q[15];
x q[15];
sdg q[17];
sdg q[18];
cx q[17], q[18];
sdg q[17];
s q[25];
s q[27];
cx q[29], q[30];
h q[31];
s q[35];
h q[39];
t q[29];
h q[20];
s q[35];
t q[30];
s q[22];
ccx q[20], q[22], q[35];
ccx q[27], q[26], q[38];
ccx q[22], q[26], q[20];
h q[23];
s q[29];
t q[22];
s q[20];
h q[34];
cx q[32], q[26];
cx q[40], q[18];
cx q[41], q[34];
cx q[42], q[7];
cx q[43], q[10];
cx q[44], q[13];
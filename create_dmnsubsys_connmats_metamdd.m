% This script extracts DMN matrices from Power atlas timeseries of the Meta-MDD
% dataset

%%% CREATE A STACK OF DMN MATRICES

clear

% List of subjects from Chao-Gan et al.
subs={'S1-1-0001', 'S1-1-0002', 'S1-1-0003', 'S1-1-0004', 'S1-1-0005', 'S1-1-0006', 'S1-1-0007', 'S1-1-0008', 'S1-1-0009', 'S1-1-0010', 'S1-1-0011', 'S1-1-0012', 'S1-1-0013', 'S1-1-0014', 'S1-1-0015', 'S1-1-0016', 'S1-1-0017', 'S1-1-0018', 'S1-1-0019', 'S1-1-0020', 'S1-1-0021', 'S1-1-0022', 'S1-1-0023', 'S1-1-0024', 'S1-1-0025', 'S1-1-0026', 'S1-1-0027', 'S1-1-0028', 'S1-1-0029', 'S1-1-0030', 'S1-1-0031', 'S1-1-0032', 'S1-1-0033', 'S1-1-0034', 'S1-1-0035', 'S1-1-0036', 'S1-1-0037', 'S1-1-0038', 'S1-1-0039', 'S1-1-0040', 'S1-1-0041', 'S1-1-0042', 'S1-1-0043', 'S1-1-0044', 'S1-1-0045', 'S1-1-0046', 'S1-1-0047', 'S1-1-0048', 'S1-1-0049', 'S1-1-0050', 'S1-1-0051', 'S1-1-0052', 'S1-1-0053', 'S1-1-0054', 'S1-1-0055', 'S1-1-0056', 'S1-1-0057', 'S1-1-0058', 'S1-1-0059', 'S1-1-0060', 'S1-1-0061', 'S1-1-0062', 'S1-1-0063', 'S1-1-0064', 'S1-1-0065', 'S1-1-0066', 'S1-1-0067', 'S1-1-0068', 'S1-1-0069', 'S1-1-0070', 'S1-1-0071', 'S1-1-0072', 'S1-1-0073', 'S1-1-0074', 'S2-1-0001', 'S2-1-0002', 'S2-1-0003', 'S2-1-0004', 'S2-1-0005', 'S2-1-0006', 'S2-1-0007', 'S2-1-0008', 'S2-1-0009', 'S2-1-0010', 'S2-1-0011', 'S2-1-0012', 'S2-1-0013', 'S2-1-0014', 'S2-1-0015', 'S2-1-0016', 'S2-1-0017', 'S2-1-0018', 'S2-1-0019', 'S2-1-0020', 'S2-1-0021', 'S2-1-0022', 'S2-1-0023', 'S2-1-0024', 'S2-1-0025', 'S2-1-0026', 'S2-1-0027', 'S2-1-0028', 'S2-1-0029', 'S2-1-0030', 'S3-1-0001', 'S3-1-0002', 'S3-1-0003', 'S3-1-0004', 'S3-1-0005', 'S3-1-0006', 'S3-1-0007', 'S3-1-0008', 'S3-1-0009', 'S3-1-0010', 'S3-1-0011', 'S3-1-0012', 'S3-1-0013', 'S3-1-0014', 'S3-1-0015', 'S3-1-0016', 'S3-1-0017', 'S3-1-0018', 'S3-1-0019', 'S3-1-0020', 'S3-1-0021', 'S3-1-0022', 'S3-1-0023', 'S3-1-0024', 'S3-1-0025', 'S3-1-0026', 'S3-1-0027', 'S5-1-0001', 'S5-1-0002', 'S5-1-0003', 'S5-1-0004', 'S5-1-0005', 'S5-1-0006', 'S5-1-0007', 'S5-1-0008', 'S5-1-0009', 'S5-1-0010', 'S5-1-0011', 'S5-1-0012', 'S5-1-0013', 'S6-1-0001', 'S6-1-0002', 'S6-1-0003', 'S6-1-0004', 'S6-1-0005', 'S6-1-0006', 'S6-1-0007', 'S6-1-0008', 'S6-1-0009', 'S6-1-0010', 'S6-1-0011', 'S6-1-0012', 'S6-1-0013', 'S6-1-0014', 'S6-1-0015', 'S7-1-0001', 'S7-1-0002', 'S7-1-0003', 'S7-1-0004', 'S7-1-0005', 'S7-1-0006', 'S7-1-0007', 'S7-1-0008', 'S7-1-0009', 'S7-1-0010', 'S7-1-0011', 'S7-1-0012', 'S7-1-0013', 'S7-1-0014', 'S7-1-0015', 'S7-1-0016', 'S7-1-0017', 'S7-1-0018', 'S7-1-0019', 'S7-1-0020', 'S7-1-0021', 'S7-1-0022', 'S7-1-0023', 'S7-1-0024', 'S7-1-0025', 'S7-1-0026', 'S7-1-0027', 'S7-1-0028', 'S7-1-0029', 'S7-1-0030', 'S7-1-0031', 'S7-1-0032', 'S7-1-0033', 'S7-1-0034', 'S7-1-0035', 'S7-1-0036', 'S7-1-0037', 'S7-1-0038', 'S8-1-0001', 'S8-1-0002', 'S8-1-0003', 'S8-1-0004', 'S8-1-0005', 'S8-1-0006', 'S8-1-0007', 'S8-1-0008', 'S8-1-0009', 'S8-1-0010', 'S8-1-0011', 'S8-1-0012', 'S8-1-0013', 'S8-1-0014', 'S8-1-0015', 'S8-1-0016', 'S8-1-0017', 'S8-1-0018', 'S8-1-0019', 'S8-1-0020', 'S8-1-0021', 'S8-1-0022', 'S8-1-0023', 'S8-1-0024', 'S8-1-0025', 'S8-1-0026', 'S8-1-0027', 'S8-1-0028', 'S8-1-0029', 'S8-1-0030', 'S8-1-0031', 'S8-1-0032', 'S8-1-0033', 'S8-1-0034', 'S8-1-0035', 'S8-1-0036', 'S8-1-0037', 'S8-1-0038', 'S8-1-0039', 'S8-1-0040', 'S8-1-0041', 'S8-1-0042', 'S8-1-0043', 'S8-1-0044', 'S8-1-0045', 'S8-1-0046', 'S8-1-0047', 'S8-1-0048', 'S8-1-0049', 'S8-1-0050', 'S8-1-0051', 'S8-1-0052', 'S8-1-0053', 'S8-1-0054', 'S8-1-0055', 'S8-1-0056', 'S8-1-0057', 'S8-1-0058', 'S8-1-0059', 'S8-1-0060', 'S8-1-0061', 'S8-1-0062', 'S8-1-0063', 'S8-1-0064', 'S8-1-0065', 'S8-1-0066', 'S8-1-0067', 'S8-1-0068', 'S8-1-0069', 'S8-1-0070', 'S8-1-0071', 'S8-1-0072', 'S8-1-0073', 'S8-1-0074', 'S8-1-0075', 'S9-1-0001', 'S9-1-0002', 'S9-1-0003', 'S9-1-0004', 'S9-1-0005', 'S9-1-0006', 'S9-1-0007', 'S9-1-0008', 'S9-1-0009', 'S9-1-0010', 'S9-1-0011', 'S9-1-0012', 'S9-1-0013', 'S9-1-0014', 'S9-1-0015', 'S9-1-0016', 'S9-1-0017', 'S9-1-0018', 'S9-1-0019', 'S9-1-0020', 'S9-1-0021', 'S9-1-0022', 'S9-1-0023', 'S9-1-0024', 'S9-1-0025', 'S9-1-0026', 'S9-1-0027', 'S9-1-0028', 'S9-1-0029', 'S9-1-0030', 'S9-1-0031', 'S9-1-0032', 'S9-1-0033', 'S9-1-0034', 'S9-1-0035', 'S9-1-0036', 'S9-1-0037', 'S9-1-0038', 'S9-1-0039', 'S9-1-0040', 'S9-1-0041', 'S9-1-0042', 'S9-1-0043', 'S9-1-0044', 'S9-1-0045', 'S9-1-0046', 'S9-1-0047', 'S9-1-0048', 'S9-1-0049', 'S9-1-0050', 'S10-1-0001', 'S10-1-0002', 'S10-1-0003', 'S10-1-0004', 'S10-1-0005', 'S10-1-0006', 'S10-1-0007', 'S10-1-0008', 'S10-1-0009', 'S10-1-0010', 'S10-1-0011', 'S10-1-0012', 'S10-1-0013', 'S10-1-0014', 'S10-1-0015', 'S10-1-0016', 'S10-1-0017', 'S10-1-0018', 'S10-1-0019', 'S10-1-0020', 'S10-1-0021', 'S10-1-0022', 'S10-1-0023', 'S10-1-0024', 'S10-1-0025', 'S10-1-0026', 'S10-1-0027', 'S10-1-0028', 'S10-1-0029', 'S10-1-0030', 'S10-1-0031', 'S10-1-0032', 'S10-1-0033', 'S10-1-0034', 'S10-1-0035', 'S10-1-0036', 'S10-1-0037', 'S10-1-0038', 'S10-1-0039', 'S10-1-0040', 'S10-1-0041', 'S10-1-0042', 'S10-1-0043', 'S10-1-0044', 'S10-1-0045', 'S10-1-0046', 'S10-1-0047', 'S10-1-0048', 'S10-1-0049', 'S10-1-0050', 'S11-1-0001', 'S11-1-0002', 'S11-1-0003', 'S11-1-0004', 'S11-1-0005', 'S11-1-0006', 'S11-1-0007', 'S11-1-0008', 'S11-1-0009', 'S11-1-0010', 'S11-1-0011', 'S11-1-0012', 'S11-1-0013', 'S11-1-0014', 'S11-1-0015', 'S11-1-0016', 'S11-1-0017', 'S11-1-0018', 'S11-1-0019', 'S11-1-0020', 'S11-1-0021', 'S11-1-0022', 'S11-1-0023', 'S11-1-0024', 'S11-1-0025', 'S11-1-0026', 'S11-1-0027', 'S11-1-0028', 'S11-1-0029', 'S11-1-0030', 'S11-1-0031', 'S11-1-0032', 'S12-1-0001', 'S12-1-0002', 'S12-1-0003', 'S12-1-0004', 'S12-1-0005', 'S12-1-0006', 'S12-1-0007', 'S12-1-0008', 'S12-1-0009', 'S12-1-0010', 'S12-1-0011', 'S12-1-0012', 'S12-1-0013', 'S12-1-0014', 'S12-1-0015', 'S12-1-0016', 'S12-1-0017', 'S12-1-0018', 'S12-1-0019', 'S12-1-0020', 'S12-1-0021', 'S12-1-0022', 'S12-1-0023', 'S12-1-0024', 'S12-1-0025', 'S12-1-0026', 'S12-1-0027', 'S12-1-0028', 'S12-1-0029', 'S12-1-0030', 'S12-1-0031', 'S12-1-0032', 'S13-1-0001', 'S13-1-0002', 'S13-1-0003', 'S13-1-0004', 'S13-1-0005', 'S13-1-0006', 'S13-1-0007', 'S13-1-0008', 'S13-1-0009', 'S13-1-0010', 'S13-1-0011', 'S13-1-0012', 'S13-1-0013', 'S13-1-0014', 'S13-1-0015', 'S13-1-0016', 'S13-1-0017', 'S13-1-0018', 'S13-1-0019', 'S13-1-0020', 'S13-1-0021', 'S13-1-0022', 'S13-1-0023', 'S13-1-0024', 'S13-1-0025', 'S14-1-0001', 'S14-1-0002', 'S14-1-0003', 'S14-1-0004', 'S14-1-0005', 'S14-1-0006', 'S14-1-0007', 'S14-1-0008', 'S14-1-0009', 'S14-1-0010', 'S14-1-0011', 'S14-1-0012', 'S14-1-0013', 'S14-1-0014', 'S14-1-0015', 'S14-1-0016', 'S14-1-0017', 'S14-1-0018', 'S14-1-0019', 'S14-1-0020', 'S14-1-0021', 'S14-1-0022', 'S14-1-0023', 'S14-1-0024', 'S14-1-0025', 'S14-1-0026', 'S14-1-0027', 'S14-1-0028', 'S14-1-0029', 'S14-1-0030', 'S14-1-0031', 'S14-1-0032', 'S14-1-0033', 'S14-1-0034', 'S14-1-0035', 'S14-1-0036', 'S14-1-0037', 'S14-1-0038', 'S14-1-0039', 'S14-1-0040', 'S14-1-0041', 'S14-1-0042', 'S14-1-0043', 'S14-1-0044', 'S14-1-0045', 'S14-1-0046', 'S14-1-0047', 'S14-1-0048', 'S14-1-0049', 'S14-1-0050', 'S14-1-0051', 'S14-1-0052', 'S14-1-0053', 'S14-1-0054', 'S14-1-0055', 'S14-1-0056', 'S14-1-0057', 'S14-1-0058', 'S14-1-0059', 'S14-1-0060', 'S14-1-0061', 'S14-1-0062', 'S14-1-0063', 'S14-1-0064', 'S15-1-0001', 'S15-1-0002', 'S15-1-0003', 'S15-1-0004', 'S15-1-0005', 'S15-1-0006', 'S15-1-0007', 'S15-1-0008', 'S15-1-0009', 'S15-1-0010', 'S15-1-0011', 'S15-1-0012', 'S15-1-0013', 'S15-1-0014', 'S15-1-0015', 'S15-1-0016', 'S15-1-0017', 'S15-1-0018', 'S15-1-0019', 'S15-1-0020', 'S15-1-0021', 'S15-1-0022', 'S15-1-0023', 'S15-1-0024', 'S15-1-0025', 'S15-1-0026', 'S15-1-0027', 'S15-1-0028', 'S15-1-0029', 'S15-1-0030', 'S15-1-0031', 'S15-1-0032', 'S15-1-0033', 'S15-1-0034', 'S15-1-0035', 'S15-1-0036', 'S15-1-0037', 'S15-1-0038', 'S15-1-0039', 'S15-1-0040', 'S15-1-0041', 'S15-1-0042', 'S15-1-0043', 'S15-1-0044', 'S15-1-0045', 'S15-1-0046', 'S15-1-0047', 'S15-1-0048', 'S15-1-0049', 'S15-1-0050', 'S16-1-0001', 'S16-1-0002', 'S16-1-0003', 'S16-1-0004', 'S16-1-0005', 'S16-1-0006', 'S16-1-0007', 'S16-1-0008', 'S16-1-0009', 'S16-1-0010', 'S16-1-0011', 'S16-1-0012', 'S16-1-0013', 'S16-1-0014', 'S16-1-0015', 'S16-1-0016', 'S16-1-0017', 'S16-1-0018', 'S16-1-0019', 'S16-1-0020', 'S16-1-0021', 'S16-1-0022', 'S16-1-0023', 'S16-1-0024', 'S16-1-0025', 'S16-1-0026', 'S16-1-0027', 'S16-1-0028', 'S16-1-0029', 'S16-1-0030', 'S16-1-0031', 'S17-1-0001', 'S17-1-0002', 'S17-1-0003', 'S17-1-0004', 'S17-1-0005', 'S17-1-0006', 'S17-1-0007', 'S17-1-0008', 'S17-1-0009', 'S17-1-0010', 'S17-1-0011', 'S17-1-0012', 'S17-1-0013', 'S17-1-0014', 'S17-1-0015', 'S17-1-0016', 'S17-1-0017', 'S17-1-0018', 'S17-1-0019', 'S17-1-0020', 'S17-1-0021', 'S17-1-0022', 'S17-1-0023', 'S17-1-0024', 'S17-1-0025', 'S17-1-0026', 'S17-1-0027', 'S17-1-0028', 'S17-1-0029', 'S17-1-0030', 'S17-1-0031', 'S17-1-0032', 'S17-1-0033', 'S17-1-0034', 'S17-1-0035', 'S17-1-0036', 'S17-1-0037', 'S17-1-0038', 'S17-1-0039', 'S17-1-0040', 'S17-1-0041', 'S17-1-0042', 'S17-1-0043', 'S17-1-0044', 'S17-1-0045', 'S17-1-0046', 'S17-1-0047', 'S18-1-0001', 'S18-1-0002', 'S18-1-0003', 'S18-1-0004', 'S18-1-0005', 'S18-1-0006', 'S18-1-0007', 'S18-1-0008', 'S18-1-0009', 'S18-1-0010', 'S18-1-0011', 'S18-1-0012', 'S18-1-0013', 'S18-1-0014', 'S18-1-0015', 'S18-1-0016', 'S18-1-0017', 'S18-1-0018', 'S18-1-0019', 'S18-1-0020', 'S18-1-0021', 'S19-1-0001', 'S19-1-0002', 'S19-1-0003', 'S19-1-0004', 'S19-1-0005', 'S19-1-0006', 'S19-1-0007', 'S19-1-0008', 'S19-1-0009', 'S19-1-0010', 'S19-1-0011', 'S19-1-0012', 'S19-1-0013', 'S19-1-0014', 'S19-1-0015', 'S19-1-0016', 'S19-1-0017', 'S19-1-0018', 'S19-1-0019', 'S19-1-0020', 'S19-1-0021', 'S19-1-0022', 'S19-1-0023', 'S19-1-0024', 'S19-1-0025', 'S19-1-0026', 'S19-1-0027', 'S19-1-0028', 'S19-1-0029', 'S19-1-0030', 'S19-1-0031', 'S19-1-0032', 'S19-1-0033', 'S19-1-0034', 'S19-1-0035', 'S19-1-0036', 'S19-1-0037', 'S19-1-0038', 'S19-1-0039', 'S19-1-0040', 'S19-1-0041', 'S19-1-0042', 'S19-1-0043', 'S19-1-0044', 'S19-1-0045', 'S19-1-0046', 'S19-1-0047', 'S19-1-0048', 'S19-1-0049', 'S19-1-0050', 'S19-1-0051', 'S20-1-0001', 'S20-1-0002', 'S20-1-0003', 'S20-1-0004', 'S20-1-0005', 'S20-1-0006', 'S20-1-0007', 'S20-1-0008', 'S20-1-0009', 'S20-1-0010', 'S20-1-0011', 'S20-1-0012', 'S20-1-0013', 'S20-1-0014', 'S20-1-0015', 'S20-1-0016', 'S20-1-0017', 'S20-1-0018', 'S20-1-0019', 'S20-1-0020', 'S20-1-0021', 'S20-1-0022', 'S20-1-0023', 'S20-1-0024', 'S20-1-0025', 'S20-1-0026', 'S20-1-0027', 'S20-1-0028', 'S20-1-0029', 'S20-1-0030', 'S20-1-0031', 'S20-1-0032', 'S20-1-0033', 'S20-1-0034', 'S20-1-0035', 'S20-1-0036', 'S20-1-0037', 'S20-1-0038', 'S20-1-0039', 'S20-1-0040', 'S20-1-0041', 'S20-1-0042', 'S20-1-0043', 'S20-1-0044', 'S20-1-0045', 'S20-1-0046', 'S20-1-0047', 'S20-1-0048', 'S20-1-0049', 'S20-1-0050', 'S20-1-0051', 'S20-1-0052', 'S20-1-0053', 'S20-1-0054', 'S20-1-0055', 'S20-1-0056', 'S20-1-0057', 'S20-1-0058', 'S20-1-0059', 'S20-1-0060', 'S20-1-0061', 'S20-1-0062', 'S20-1-0063', 'S20-1-0064', 'S20-1-0065', 'S20-1-0066', 'S20-1-0067', 'S20-1-0068', 'S20-1-0069', 'S20-1-0070', 'S20-1-0071', 'S20-1-0072', 'S20-1-0073', 'S20-1-0074', 'S20-1-0075', 'S20-1-0076', 'S20-1-0077', 'S20-1-0078', 'S20-1-0079', 'S20-1-0080', 'S20-1-0081', 'S20-1-0082', 'S20-1-0083', 'S20-1-0084', 'S20-1-0085', 'S20-1-0086', 'S20-1-0087', 'S20-1-0088', 'S20-1-0089', 'S20-1-0090', 'S20-1-0091', 'S20-1-0092', 'S20-1-0093', 'S20-1-0094', 'S20-1-0095', 'S20-1-0096', 'S20-1-0097', 'S20-1-0098', 'S20-1-0099', 'S20-1-0100', 'S20-1-0101', 'S20-1-0102', 'S20-1-0103', 'S20-1-0104', 'S20-1-0105', 'S20-1-0106', 'S20-1-0107', 'S20-1-0108', 'S20-1-0109', 'S20-1-0110', 'S20-1-0111', 'S20-1-0112', 'S20-1-0113', 'S20-1-0114', 'S20-1-0115', 'S20-1-0116', 'S20-1-0117', 'S20-1-0118', 'S20-1-0119', 'S20-1-0120', 'S20-1-0121', 'S20-1-0122', 'S20-1-0123', 'S20-1-0124', 'S20-1-0125', 'S20-1-0126', 'S20-1-0127', 'S20-1-0128', 'S20-1-0129', 'S20-1-0130', 'S20-1-0131', 'S20-1-0132', 'S20-1-0133', 'S20-1-0134', 'S20-1-0135', 'S20-1-0136', 'S20-1-0137', 'S20-1-0138', 'S20-1-0139', 'S20-1-0140', 'S20-1-0141', 'S20-1-0142', 'S20-1-0143', 'S20-1-0144', 'S20-1-0145', 'S20-1-0146', 'S20-1-0147', 'S20-1-0148', 'S20-1-0149', 'S20-1-0150', 'S20-1-0151', 'S20-1-0152', 'S20-1-0153', 'S20-1-0154', 'S20-1-0155', 'S20-1-0156', 'S20-1-0157', 'S20-1-0158', 'S20-1-0159', 'S20-1-0160', 'S20-1-0161', 'S20-1-0162', 'S20-1-0163', 'S20-1-0164', 'S20-1-0165', 'S20-1-0166', 'S20-1-0167', 'S20-1-0168', 'S20-1-0169', 'S20-1-0170', 'S20-1-0171', 'S20-1-0172', 'S20-1-0173', 'S20-1-0174', 'S20-1-0175', 'S20-1-0176', 'S20-1-0177', 'S20-1-0178', 'S20-1-0179', 'S20-1-0180', 'S20-1-0181', 'S20-1-0182', 'S20-1-0183', 'S20-1-0184', 'S20-1-0185', 'S20-1-0186', 'S20-1-0187', 'S20-1-0188', 'S20-1-0189', 'S20-1-0190', 'S20-1-0191', 'S20-1-0192', 'S20-1-0193', 'S20-1-0194', 'S20-1-0195', 'S20-1-0196', 'S20-1-0197', 'S20-1-0198', 'S20-1-0199', 'S20-1-0200', 'S20-1-0201', 'S20-1-0202', 'S20-1-0203', 'S20-1-0204', 'S20-1-0205', 'S20-1-0206', 'S20-1-0207', 'S20-1-0208', 'S20-1-0209', 'S20-1-0210', 'S20-1-0211', 'S20-1-0212', 'S20-1-0213', 'S20-1-0214', 'S20-1-0215', 'S20-1-0216', 'S20-1-0217', 'S20-1-0218', 'S20-1-0219', 'S20-1-0220', 'S20-1-0221', 'S20-1-0222', 'S20-1-0223', 'S20-1-0224', 'S20-1-0225', 'S20-1-0226', 'S20-1-0227', 'S20-1-0228', 'S20-1-0229', 'S20-1-0230', 'S20-1-0231', 'S20-1-0232', 'S20-1-0233', 'S20-1-0234', 'S20-1-0235', 'S20-1-0236', 'S20-1-0237', 'S20-1-0238', 'S20-1-0239', 'S20-1-0240', 'S20-1-0241', 'S20-1-0242', 'S20-1-0243', 'S20-1-0244', 'S20-1-0245', 'S20-1-0246', 'S20-1-0247', 'S20-1-0248', 'S20-1-0249', 'S20-1-0250', 'S20-1-0251', 'S20-1-0252', 'S20-1-0253', 'S20-1-0254', 'S20-1-0255', 'S20-1-0256', 'S20-1-0257', 'S20-1-0258', 'S20-1-0259', 'S20-1-0260', 'S20-1-0261', 'S20-1-0262', 'S20-1-0263', 'S20-1-0264', 'S20-1-0265', 'S20-1-0266', 'S20-1-0267', 'S20-1-0268', 'S20-1-0269', 'S20-1-0270', 'S20-1-0271', 'S20-1-0272', 'S20-1-0273', 'S20-1-0274', 'S20-1-0275', 'S20-1-0276', 'S20-1-0277', 'S20-1-0278', 'S20-1-0279', 'S20-1-0280', 'S20-1-0281', 'S20-1-0282', 'S21-1-0001', 'S21-1-0002', 'S21-1-0003', 'S21-1-0004', 'S21-1-0005', 'S21-1-0006', 'S21-1-0007', 'S21-1-0008', 'S21-1-0009', 'S21-1-0010', 'S21-1-0011', 'S21-1-0012', 'S21-1-0013', 'S21-1-0014', 'S21-1-0015', 'S21-1-0016', 'S21-1-0017', 'S21-1-0018', 'S21-1-0019', 'S21-1-0020', 'S21-1-0021', 'S21-1-0022', 'S21-1-0023', 'S21-1-0024', 'S21-1-0025', 'S21-1-0026', 'S21-1-0027', 'S21-1-0028', 'S21-1-0029', 'S21-1-0030', 'S21-1-0031', 'S21-1-0032', 'S21-1-0033', 'S21-1-0034', 'S21-1-0035', 'S21-1-0036', 'S21-1-0037', 'S21-1-0038', 'S21-1-0039', 'S21-1-0040', 'S21-1-0041', 'S21-1-0042', 'S21-1-0043', 'S21-1-0044', 'S21-1-0045', 'S21-1-0046', 'S21-1-0047', 'S21-1-0048', 'S21-1-0049', 'S21-1-0050', 'S21-1-0051', 'S21-1-0052', 'S21-1-0053', 'S21-1-0054', 'S21-1-0055', 'S21-1-0056', 'S21-1-0057', 'S21-1-0058', 'S21-1-0059', 'S21-1-0060', 'S21-1-0061', 'S21-1-0062', 'S21-1-0063', 'S21-1-0064', 'S21-1-0065', 'S21-1-0066', 'S21-1-0067', 'S21-1-0068', 'S21-1-0069', 'S21-1-0070', 'S21-1-0071', 'S21-1-0072', 'S21-1-0073', 'S21-1-0074', 'S21-1-0075', 'S21-1-0076', 'S21-1-0077', 'S21-1-0078', 'S21-1-0079', 'S21-1-0080', 'S21-1-0081', 'S21-1-0082', 'S21-1-0083', 'S21-1-0084', 'S21-1-0085', 'S21-1-0086', 'S22-1-0001', 'S22-1-0002', 'S22-1-0003', 'S22-1-0004', 'S22-1-0005', 'S22-1-0006', 'S22-1-0007', 'S22-1-0008', 'S22-1-0009', 'S22-1-0010', 'S22-1-0011', 'S22-1-0012', 'S22-1-0013', 'S22-1-0014', 'S22-1-0015', 'S22-1-0016', 'S22-1-0017', 'S22-1-0018', 'S22-1-0019', 'S22-1-0020', 'S22-1-0021', 'S22-1-0022', 'S22-1-0023', 'S22-1-0024', 'S22-1-0025', 'S22-1-0026', 'S22-1-0027', 'S22-1-0028', 'S22-1-0029', 'S22-1-0030', 'S23-1-0001', 'S23-1-0002', 'S23-1-0003', 'S23-1-0004', 'S23-1-0005', 'S23-1-0006', 'S23-1-0007', 'S23-1-0008', 'S23-1-0009', 'S23-1-0010', 'S23-1-0011', 'S23-1-0012', 'S23-1-0013', 'S23-1-0014', 'S23-1-0015', 'S23-1-0016', 'S23-1-0017', 'S23-1-0018', 'S23-1-0019', 'S23-1-0020', 'S23-1-0021', 'S23-1-0022', 'S23-1-0023', 'S23-1-0024', 'S23-1-0025', 'S23-1-0026', 'S23-1-0027', 'S23-1-0028', 'S23-1-0029', 'S23-1-0030', 'S23-1-0031', 'S23-1-0032', 'S24-1-0001', 'S24-1-0002', 'S24-1-0003', 'S24-1-0004', 'S24-1-0005', 'S24-1-0006', 'S24-1-0007', 'S24-1-0008', 'S24-1-0009', 'S24-1-0010', 'S24-1-0011', 'S24-1-0012', 'S24-1-0013', 'S24-1-0014', 'S24-1-0015', 'S24-1-0016', 'S24-1-0017', 'S24-1-0018', 'S24-1-0019', 'S24-1-0020', 'S24-1-0021', 'S24-1-0022', 'S24-1-0023', 'S24-1-0024', 'S24-1-0025', 'S24-1-0026', 'S24-1-0027', 'S24-1-0028', 'S24-1-0029', 'S24-1-0030', 'S24-1-0031', 'S24-1-0032', 'S25-1-0001', 'S25-1-0002', 'S25-1-0003', 'S25-1-0004', 'S25-1-0005', 'S25-1-0006', 'S25-1-0007', 'S25-1-0008', 'S25-1-0009', 'S25-1-0010', 'S25-1-0011', 'S25-1-0012', 'S25-1-0013', 'S25-1-0014', 'S25-1-0015', 'S25-1-0016', 'S25-1-0017', 'S25-1-0018', 'S25-1-0019', 'S25-1-0020', 'S25-1-0021', 'S25-1-0022', 'S25-1-0023', 'S25-1-0024', 'S25-1-0025', 'S25-1-0026', 'S25-1-0027', 'S25-1-0028', 'S25-1-0029', 'S25-1-0030', 'S25-1-0031', 'S25-1-0032', 'S25-1-0033', 'S25-1-0034', 'S25-1-0035', 'S25-1-0036', 'S25-1-0037', 'S25-1-0038', 'S25-1-0039', 'S25-1-0040', 'S25-1-0041', 'S25-1-0042', 'S25-1-0043', 'S25-1-0044', 'S25-1-0045', 'S25-1-0046', 'S25-1-0047', 'S25-1-0048', 'S25-1-0049', 'S25-1-0050', 'S25-1-0051', 'S25-1-0052', 'S25-1-0053', 'S25-1-0054', 'S25-1-0055', 'S25-1-0056', 'S25-1-0057', 'S25-1-0058', 'S25-1-0059', 'S25-1-0060', 'S25-1-0061', 'S25-1-0062', 'S25-1-0063', 'S25-1-0064', 'S25-1-0065', 'S25-1-0066', 'S25-1-0067', 'S25-1-0068', 'S25-1-0069', 'S25-1-0070', 'S25-1-0071', 'S25-1-0072', 'S25-1-0073', 'S25-1-0074', 'S25-1-0075', 'S25-1-0076', 'S25-1-0077', 'S25-1-0078', 'S25-1-0079', 'S25-1-0080', 'S25-1-0081', 'S25-1-0082', 'S25-1-0083', 'S25-1-0084', 'S25-1-0085', 'S25-1-0086', 'S25-1-0087', 'S25-1-0088', 'S25-1-0089', 'S1-2-0001', 'S1-2-0002', 'S1-2-0003', 'S1-2-0004', 'S1-2-0005', 'S1-2-0006', 'S1-2-0007', 'S1-2-0008', 'S1-2-0009', 'S1-2-0010', 'S1-2-0011', 'S1-2-0012', 'S1-2-0013', 'S1-2-0014', 'S1-2-0015', 'S1-2-0016', 'S1-2-0017', 'S1-2-0018', 'S1-2-0019', 'S1-2-0020', 'S1-2-0021', 'S1-2-0022', 'S1-2-0023', 'S1-2-0024', 'S1-2-0025', 'S1-2-0026', 'S1-2-0027', 'S1-2-0028', 'S1-2-0029', 'S1-2-0030', 'S1-2-0031', 'S1-2-0032', 'S1-2-0033', 'S1-2-0034', 'S1-2-0035', 'S1-2-0036', 'S1-2-0037', 'S1-2-0038', 'S1-2-0039', 'S1-2-0040', 'S1-2-0041', 'S1-2-0042', 'S1-2-0043', 'S1-2-0044', 'S1-2-0045', 'S1-2-0046', 'S1-2-0047', 'S1-2-0048', 'S1-2-0049', 'S1-2-0050', 'S1-2-0051', 'S1-2-0052', 'S1-2-0053', 'S1-2-0054', 'S1-2-0055', 'S1-2-0056', 'S1-2-0057', 'S1-2-0058', 'S1-2-0059', 'S1-2-0060', 'S1-2-0061', 'S1-2-0062', 'S1-2-0063', 'S1-2-0064', 'S1-2-0065', 'S1-2-0066', 'S1-2-0067', 'S1-2-0068', 'S1-2-0069', 'S1-2-0070', 'S1-2-0071', 'S1-2-0072', 'S1-2-0073', 'S1-2-0074', 'S2-2-0001', 'S2-2-0002', 'S2-2-0003', 'S2-2-0004', 'S2-2-0005', 'S2-2-0006', 'S2-2-0007', 'S2-2-0008', 'S2-2-0009', 'S2-2-0010', 'S2-2-0011', 'S2-2-0012', 'S2-2-0013', 'S2-2-0014', 'S2-2-0015', 'S2-2-0016', 'S2-2-0017', 'S2-2-0018', 'S2-2-0019', 'S2-2-0020', 'S2-2-0021', 'S2-2-0022', 'S2-2-0023', 'S2-2-0024', 'S2-2-0025', 'S2-2-0026', 'S2-2-0027', 'S2-2-0028', 'S2-2-0029', 'S2-2-0030', 'S3-2-0001', 'S3-2-0002', 'S3-2-0003', 'S3-2-0004', 'S3-2-0005', 'S3-2-0006', 'S3-2-0007', 'S3-2-0008', 'S3-2-0009', 'S3-2-0010', 'S3-2-0011', 'S3-2-0012', 'S3-2-0013', 'S3-2-0014', 'S3-2-0015', 'S3-2-0016', 'S3-2-0017', 'S3-2-0018', 'S3-2-0019', 'S3-2-0020', 'S3-2-0021', 'S3-2-0022', 'S3-2-0023', 'S3-2-0024', 'S3-2-0025', 'S3-2-0026', 'S3-2-0027', 'S3-2-0028', 'S3-2-0029', 'S3-2-0030', 'S3-2-0031', 'S3-2-0032', 'S3-2-0033', 'S3-2-0034', 'S3-2-0035', 'S3-2-0036', 'S3-2-0037', 'S5-2-0001', 'S5-2-0002', 'S5-2-0003', 'S5-2-0004', 'S5-2-0005', 'S5-2-0006', 'S5-2-0007', 'S5-2-0008', 'S5-2-0009', 'S5-2-0010', 'S5-2-0011', 'S6-2-0001', 'S6-2-0002', 'S6-2-0003', 'S6-2-0004', 'S6-2-0005', 'S6-2-0006', 'S6-2-0007', 'S6-2-0008', 'S6-2-0009', 'S6-2-0010', 'S6-2-0011', 'S6-2-0012', 'S6-2-0013', 'S6-2-0014', 'S6-2-0015', 'S7-2-0001', 'S7-2-0002', 'S7-2-0003', 'S7-2-0004', 'S7-2-0005', 'S7-2-0006', 'S7-2-0007', 'S7-2-0008', 'S7-2-0009', 'S7-2-0010', 'S7-2-0011', 'S7-2-0012', 'S7-2-0013', 'S7-2-0014', 'S7-2-0015', 'S7-2-0016', 'S7-2-0017', 'S7-2-0018', 'S7-2-0019', 'S7-2-0020', 'S7-2-0021', 'S7-2-0022', 'S7-2-0023', 'S7-2-0024', 'S7-2-0025', 'S7-2-0026', 'S7-2-0027', 'S7-2-0028', 'S7-2-0029', 'S7-2-0030', 'S7-2-0031', 'S7-2-0032', 'S7-2-0033', 'S7-2-0034', 'S7-2-0035', 'S7-2-0036', 'S7-2-0037', 'S7-2-0038', 'S7-2-0039', 'S7-2-0040', 'S7-2-0041', 'S7-2-0042', 'S7-2-0043', 'S7-2-0044', 'S7-2-0045', 'S7-2-0046', 'S7-2-0047', 'S7-2-0048', 'S7-2-0049', 'S8-2-0001', 'S8-2-0002', 'S8-2-0003', 'S8-2-0004', 'S8-2-0005', 'S8-2-0006', 'S8-2-0007', 'S8-2-0008', 'S8-2-0009', 'S8-2-0010', 'S8-2-0011', 'S8-2-0012', 'S8-2-0013', 'S8-2-0014', 'S8-2-0015', 'S8-2-0016', 'S8-2-0017', 'S8-2-0018', 'S8-2-0019', 'S8-2-0020', 'S8-2-0021', 'S8-2-0022', 'S8-2-0023', 'S8-2-0024', 'S8-2-0025', 'S8-2-0026', 'S8-2-0027', 'S8-2-0028', 'S8-2-0029', 'S8-2-0030', 'S8-2-0031', 'S8-2-0032', 'S8-2-0033', 'S8-2-0034', 'S8-2-0035', 'S8-2-0036', 'S8-2-0037', 'S8-2-0038', 'S8-2-0039', 'S8-2-0040', 'S8-2-0041', 'S8-2-0042', 'S8-2-0043', 'S8-2-0044', 'S8-2-0045', 'S8-2-0046', 'S8-2-0047', 'S8-2-0048', 'S8-2-0049', 'S8-2-0050', 'S8-2-0051', 'S8-2-0052', 'S8-2-0053', 'S8-2-0054', 'S8-2-0055', 'S8-2-0056', 'S8-2-0057', 'S8-2-0058', 'S8-2-0059', 'S8-2-0060', 'S8-2-0061', 'S8-2-0062', 'S8-2-0063', 'S8-2-0064', 'S8-2-0065', 'S8-2-0066', 'S8-2-0067', 'S8-2-0068', 'S8-2-0069', 'S8-2-0070', 'S8-2-0071', 'S8-2-0072', 'S8-2-0073', 'S8-2-0074', 'S8-2-0075', 'S9-2-0001', 'S9-2-0002', 'S9-2-0003', 'S9-2-0004', 'S9-2-0005', 'S9-2-0006', 'S9-2-0007', 'S9-2-0008', 'S9-2-0009', 'S9-2-0010', 'S9-2-0011', 'S9-2-0012', 'S9-2-0013', 'S9-2-0014', 'S9-2-0015', 'S9-2-0016', 'S9-2-0017', 'S9-2-0018', 'S9-2-0019', 'S9-2-0020', 'S9-2-0021', 'S9-2-0022', 'S9-2-0023', 'S9-2-0024', 'S9-2-0025', 'S9-2-0026', 'S9-2-0027', 'S9-2-0028', 'S9-2-0029', 'S9-2-0030', 'S9-2-0031', 'S9-2-0032', 'S9-2-0033', 'S9-2-0034', 'S9-2-0035', 'S9-2-0036', 'S9-2-0037', 'S9-2-0038', 'S9-2-0039', 'S9-2-0040', 'S9-2-0041', 'S9-2-0042', 'S9-2-0043', 'S9-2-0044', 'S9-2-0045', 'S9-2-0046', 'S9-2-0047', 'S9-2-0048', 'S9-2-0049', 'S9-2-0050', 'S10-2-0001', 'S10-2-0002', 'S10-2-0003', 'S10-2-0004', 'S10-2-0005', 'S10-2-0006', 'S10-2-0007', 'S10-2-0008', 'S10-2-0009', 'S10-2-0010', 'S10-2-0011', 'S10-2-0012', 'S10-2-0013', 'S10-2-0014', 'S10-2-0015', 'S10-2-0016', 'S10-2-0017', 'S10-2-0018', 'S10-2-0019', 'S10-2-0020', 'S10-2-0021', 'S10-2-0022', 'S10-2-0023', 'S10-2-0024', 'S10-2-0025', 'S10-2-0026', 'S10-2-0027', 'S10-2-0028', 'S10-2-0029', 'S10-2-0030', 'S10-2-0031', 'S10-2-0032', 'S10-2-0033', 'S11-2-0001', 'S11-2-0002', 'S11-2-0003', 'S11-2-0004', 'S11-2-0005', 'S11-2-0006', 'S11-2-0007', 'S11-2-0008', 'S11-2-0009', 'S11-2-0010', 'S11-2-0011', 'S11-2-0012', 'S11-2-0013', 'S11-2-0014', 'S11-2-0015', 'S11-2-0016', 'S11-2-0017', 'S11-2-0018', 'S11-2-0019', 'S11-2-0020', 'S11-2-0021', 'S11-2-0022', 'S11-2-0023', 'S11-2-0024', 'S11-2-0025', 'S11-2-0026', 'S11-2-0027', 'S11-2-0028', 'S11-2-0029', 'S12-2-0001', 'S12-2-0002', 'S12-2-0003', 'S12-2-0004', 'S12-2-0005', 'S12-2-0006', 'S13-2-0001', 'S13-2-0002', 'S13-2-0003', 'S13-2-0004', 'S13-2-0005', 'S13-2-0006', 'S13-2-0007', 'S13-2-0008', 'S13-2-0009', 'S13-2-0010', 'S13-2-0011', 'S13-2-0012', 'S13-2-0013', 'S13-2-0014', 'S13-2-0015', 'S13-2-0016', 'S13-2-0017', 'S14-2-0001', 'S14-2-0002', 'S14-2-0003', 'S14-2-0004', 'S14-2-0005', 'S14-2-0006', 'S14-2-0007', 'S14-2-0008', 'S14-2-0009', 'S14-2-0010', 'S14-2-0011', 'S14-2-0012', 'S14-2-0013', 'S14-2-0014', 'S14-2-0015', 'S14-2-0016', 'S14-2-0017', 'S14-2-0018', 'S14-2-0019', 'S14-2-0020', 'S14-2-0021', 'S14-2-0022', 'S14-2-0023', 'S14-2-0024', 'S14-2-0025', 'S14-2-0026', 'S14-2-0027', 'S14-2-0028', 'S14-2-0029', 'S14-2-0030', 'S14-2-0031', 'S14-2-0032', 'S15-2-0001', 'S15-2-0002', 'S15-2-0003', 'S15-2-0004', 'S15-2-0005', 'S15-2-0006', 'S15-2-0007', 'S15-2-0008', 'S15-2-0009', 'S15-2-0010', 'S15-2-0011', 'S15-2-0012', 'S15-2-0013', 'S15-2-0014', 'S15-2-0015', 'S15-2-0016', 'S15-2-0017', 'S15-2-0018', 'S15-2-0019', 'S15-2-0020', 'S15-2-0021', 'S15-2-0022', 'S15-2-0023', 'S15-2-0024', 'S15-2-0025', 'S15-2-0026', 'S15-2-0027', 'S15-2-0028', 'S15-2-0029', 'S15-2-0030', 'S15-2-0031', 'S15-2-0032', 'S15-2-0033', 'S15-2-0034', 'S15-2-0035', 'S15-2-0036', 'S15-2-0037', 'S15-2-0038', 'S15-2-0039', 'S15-2-0040', 'S15-2-0041', 'S15-2-0042', 'S15-2-0043', 'S15-2-0044', 'S15-2-0045', 'S15-2-0046', 'S15-2-0047', 'S15-2-0048', 'S15-2-0049', 'S15-2-0050', 'S16-2-0001', 'S16-2-0002', 'S16-2-0003', 'S16-2-0004', 'S16-2-0005', 'S16-2-0006', 'S16-2-0007', 'S16-2-0008', 'S16-2-0009', 'S16-2-0010', 'S16-2-0011', 'S16-2-0012', 'S16-2-0013', 'S16-2-0014', 'S16-2-0015', 'S16-2-0016', 'S16-2-0017', 'S16-2-0018', 'S16-2-0019', 'S16-2-0020', 'S16-2-0021', 'S16-2-0022', 'S16-2-0023', 'S16-2-0024', 'S16-2-0025', 'S16-2-0026', 'S16-2-0027', 'S16-2-0028', 'S16-2-0029', 'S16-2-0030', 'S16-2-0031', 'S17-2-0001', 'S17-2-0002', 'S17-2-0003', 'S17-2-0004', 'S17-2-0005', 'S17-2-0006', 'S17-2-0007', 'S17-2-0008', 'S17-2-0009', 'S17-2-0010', 'S17-2-0011', 'S17-2-0012', 'S17-2-0013', 'S17-2-0014', 'S17-2-0015', 'S17-2-0016', 'S17-2-0017', 'S17-2-0018', 'S17-2-0019', 'S17-2-0020', 'S17-2-0021', 'S17-2-0022', 'S17-2-0023', 'S17-2-0024', 'S17-2-0025', 'S17-2-0026', 'S17-2-0027', 'S17-2-0028', 'S17-2-0029', 'S17-2-0030', 'S17-2-0031', 'S17-2-0032', 'S17-2-0033', 'S17-2-0034', 'S17-2-0035', 'S17-2-0036', 'S17-2-0037', 'S17-2-0038', 'S17-2-0039', 'S17-2-0040', 'S17-2-0041', 'S17-2-0042', 'S17-2-0043', 'S17-2-0044', 'S18-2-0001', 'S18-2-0002', 'S18-2-0003', 'S18-2-0004', 'S18-2-0005', 'S18-2-0006', 'S18-2-0007', 'S18-2-0008', 'S18-2-0009', 'S18-2-0010', 'S18-2-0011', 'S18-2-0012', 'S18-2-0013', 'S18-2-0014', 'S18-2-0015', 'S18-2-0016', 'S18-2-0017', 'S18-2-0018', 'S18-2-0019', 'S18-2-0020', 'S19-2-0001', 'S19-2-0002', 'S19-2-0003', 'S19-2-0004', 'S19-2-0005', 'S19-2-0006', 'S19-2-0007', 'S19-2-0008', 'S19-2-0009', 'S19-2-0010', 'S19-2-0011', 'S19-2-0012', 'S19-2-0013', 'S19-2-0014', 'S19-2-0015', 'S19-2-0016', 'S19-2-0017', 'S19-2-0018', 'S19-2-0019', 'S19-2-0020', 'S19-2-0021', 'S19-2-0022', 'S19-2-0023', 'S19-2-0024', 'S19-2-0025', 'S19-2-0026', 'S19-2-0027', 'S19-2-0028', 'S19-2-0029', 'S19-2-0030', 'S19-2-0031', 'S19-2-0032', 'S19-2-0033', 'S19-2-0034', 'S19-2-0035', 'S19-2-0036', 'S20-2-0001', 'S20-2-0002', 'S20-2-0003', 'S20-2-0004', 'S20-2-0005', 'S20-2-0006', 'S20-2-0007', 'S20-2-0008', 'S20-2-0009', 'S20-2-0010', 'S20-2-0011', 'S20-2-0012', 'S20-2-0013', 'S20-2-0014', 'S20-2-0015', 'S20-2-0016', 'S20-2-0017', 'S20-2-0018', 'S20-2-0019', 'S20-2-0020', 'S20-2-0021', 'S20-2-0022', 'S20-2-0023', 'S20-2-0024', 'S20-2-0025', 'S20-2-0026', 'S20-2-0027', 'S20-2-0028', 'S20-2-0029', 'S20-2-0030', 'S20-2-0031', 'S20-2-0032', 'S20-2-0033', 'S20-2-0034', 'S20-2-0035', 'S20-2-0036', 'S20-2-0037', 'S20-2-0038', 'S20-2-0039', 'S20-2-0040', 'S20-2-0041', 'S20-2-0042', 'S20-2-0043', 'S20-2-0044', 'S20-2-0045', 'S20-2-0046', 'S20-2-0047', 'S20-2-0048', 'S20-2-0049', 'S20-2-0050', 'S20-2-0051', 'S20-2-0052', 'S20-2-0053', 'S20-2-0054', 'S20-2-0055', 'S20-2-0056', 'S20-2-0057', 'S20-2-0058', 'S20-2-0059', 'S20-2-0060', 'S20-2-0061', 'S20-2-0062', 'S20-2-0063', 'S20-2-0064', 'S20-2-0065', 'S20-2-0066', 'S20-2-0067', 'S20-2-0068', 'S20-2-0069', 'S20-2-0070', 'S20-2-0071', 'S20-2-0072', 'S20-2-0073', 'S20-2-0074', 'S20-2-0075', 'S20-2-0076', 'S20-2-0077', 'S20-2-0078', 'S20-2-0079', 'S20-2-0080', 'S20-2-0081', 'S20-2-0082', 'S20-2-0083', 'S20-2-0084', 'S20-2-0085', 'S20-2-0086', 'S20-2-0087', 'S20-2-0088', 'S20-2-0089', 'S20-2-0090', 'S20-2-0091', 'S20-2-0092', 'S20-2-0093', 'S20-2-0094', 'S20-2-0095', 'S20-2-0096', 'S20-2-0097', 'S20-2-0098', 'S20-2-0099', 'S20-2-0100', 'S20-2-0101', 'S20-2-0102', 'S20-2-0103', 'S20-2-0104', 'S20-2-0105', 'S20-2-0106', 'S20-2-0107', 'S20-2-0108', 'S20-2-0109', 'S20-2-0110', 'S20-2-0111', 'S20-2-0112', 'S20-2-0113', 'S20-2-0114', 'S20-2-0115', 'S20-2-0116', 'S20-2-0117', 'S20-2-0118', 'S20-2-0119', 'S20-2-0120', 'S20-2-0121', 'S20-2-0122', 'S20-2-0123', 'S20-2-0124', 'S20-2-0125', 'S20-2-0126', 'S20-2-0127', 'S20-2-0128', 'S20-2-0129', 'S20-2-0130', 'S20-2-0131', 'S20-2-0132', 'S20-2-0133', 'S20-2-0134', 'S20-2-0135', 'S20-2-0136', 'S20-2-0137', 'S20-2-0138', 'S20-2-0139', 'S20-2-0140', 'S20-2-0141', 'S20-2-0142', 'S20-2-0143', 'S20-2-0144', 'S20-2-0145', 'S20-2-0146', 'S20-2-0147', 'S20-2-0148', 'S20-2-0149', 'S20-2-0150', 'S20-2-0151', 'S20-2-0152', 'S20-2-0153', 'S20-2-0154', 'S20-2-0155', 'S20-2-0156', 'S20-2-0157', 'S20-2-0158', 'S20-2-0159', 'S20-2-0160', 'S20-2-0161', 'S20-2-0162', 'S20-2-0163', 'S20-2-0164', 'S20-2-0165', 'S20-2-0166', 'S20-2-0167', 'S20-2-0168', 'S20-2-0169', 'S20-2-0170', 'S20-2-0171', 'S20-2-0172', 'S20-2-0173', 'S20-2-0174', 'S20-2-0175', 'S20-2-0176', 'S20-2-0177', 'S20-2-0178', 'S20-2-0179', 'S20-2-0180', 'S20-2-0181', 'S20-2-0182', 'S20-2-0183', 'S20-2-0184', 'S20-2-0185', 'S20-2-0186', 'S20-2-0187', 'S20-2-0188', 'S20-2-0189', 'S20-2-0190', 'S20-2-0191', 'S20-2-0192', 'S20-2-0193', 'S20-2-0194', 'S20-2-0195', 'S20-2-0196', 'S20-2-0197', 'S20-2-0198', 'S20-2-0199', 'S20-2-0200', 'S20-2-0201', 'S20-2-0202', 'S20-2-0203', 'S20-2-0204', 'S20-2-0205', 'S20-2-0206', 'S20-2-0207', 'S20-2-0208', 'S20-2-0209', 'S20-2-0210', 'S20-2-0211', 'S20-2-0212', 'S20-2-0213', 'S20-2-0214', 'S20-2-0215', 'S20-2-0216', 'S20-2-0217', 'S20-2-0218', 'S20-2-0219', 'S20-2-0220', 'S20-2-0221', 'S20-2-0222', 'S20-2-0223', 'S20-2-0224', 'S20-2-0225', 'S20-2-0226', 'S20-2-0227', 'S20-2-0228', 'S20-2-0229', 'S20-2-0230', 'S20-2-0231', 'S20-2-0232', 'S20-2-0233', 'S20-2-0234', 'S20-2-0235', 'S20-2-0236', 'S20-2-0237', 'S20-2-0238', 'S20-2-0239', 'S20-2-0240', 'S20-2-0241', 'S20-2-0242', 'S20-2-0243', 'S20-2-0244', 'S20-2-0245', 'S20-2-0246', 'S20-2-0247', 'S20-2-0248', 'S20-2-0249', 'S20-2-0250', 'S20-2-0251', 'S21-2-0001', 'S21-2-0002', 'S21-2-0003', 'S21-2-0004', 'S21-2-0005', 'S21-2-0006', 'S21-2-0007', 'S21-2-0008', 'S21-2-0009', 'S21-2-0010', 'S21-2-0011', 'S21-2-0012', 'S21-2-0013', 'S21-2-0014', 'S21-2-0015', 'S21-2-0016', 'S21-2-0017', 'S21-2-0018', 'S21-2-0019', 'S21-2-0020', 'S21-2-0021', 'S21-2-0022', 'S21-2-0023', 'S21-2-0024', 'S21-2-0025', 'S21-2-0026', 'S21-2-0027', 'S21-2-0028', 'S21-2-0029', 'S21-2-0030', 'S21-2-0031', 'S21-2-0032', 'S21-2-0033', 'S21-2-0034', 'S21-2-0035', 'S21-2-0036', 'S21-2-0037', 'S21-2-0038', 'S21-2-0039', 'S21-2-0040', 'S21-2-0041', 'S21-2-0042', 'S21-2-0043', 'S21-2-0044', 'S21-2-0045', 'S21-2-0046', 'S21-2-0047', 'S21-2-0048', 'S21-2-0049', 'S21-2-0050', 'S21-2-0051', 'S21-2-0052', 'S21-2-0053', 'S21-2-0054', 'S21-2-0055', 'S21-2-0056', 'S21-2-0057', 'S21-2-0058', 'S21-2-0059', 'S21-2-0060', 'S21-2-0061', 'S21-2-0062', 'S21-2-0063', 'S21-2-0064', 'S21-2-0065', 'S21-2-0066', 'S21-2-0067', 'S21-2-0068', 'S21-2-0069', 'S21-2-0070', 'S22-2-0001', 'S22-2-0002', 'S22-2-0003', 'S22-2-0004', 'S22-2-0005', 'S22-2-0006', 'S22-2-0007', 'S22-2-0008', 'S22-2-0009', 'S22-2-0010', 'S22-2-0011', 'S22-2-0012', 'S22-2-0013', 'S22-2-0014', 'S22-2-0015', 'S22-2-0016', 'S22-2-0017', 'S22-2-0018', 'S22-2-0019', 'S22-2-0020', 'S23-2-0001', 'S23-2-0002', 'S23-2-0003', 'S23-2-0004', 'S23-2-0005', 'S23-2-0006', 'S23-2-0007', 'S23-2-0008', 'S23-2-0009', 'S23-2-0010', 'S23-2-0011', 'S23-2-0012', 'S23-2-0013', 'S23-2-0014', 'S23-2-0015', 'S23-2-0016', 'S23-2-0017', 'S23-2-0018', 'S23-2-0019', 'S23-2-0020', 'S23-2-0021', 'S23-2-0022', 'S23-2-0023', 'S23-2-0024', 'S23-2-0025', 'S23-2-0026', 'S23-2-0027', 'S23-2-0028', 'S23-2-0029', 'S23-2-0030', 'S24-2-0001', 'S24-2-0002', 'S24-2-0003', 'S24-2-0004', 'S24-2-0005', 'S24-2-0006', 'S24-2-0007', 'S24-2-0008', 'S24-2-0009', 'S24-2-0010', 'S24-2-0011', 'S24-2-0012', 'S24-2-0013', 'S24-2-0014', 'S24-2-0015', 'S24-2-0016', 'S24-2-0017', 'S24-2-0018', 'S24-2-0019', 'S24-2-0020', 'S24-2-0021', 'S24-2-0022', 'S24-2-0023', 'S24-2-0024', 'S24-2-0025', 'S24-2-0026', 'S24-2-0027', 'S24-2-0028', 'S24-2-0029', 'S24-2-0030', 'S24-2-0031', 'S25-2-0001', 'S25-2-0002', 'S25-2-0003', 'S25-2-0004', 'S25-2-0005', 'S25-2-0006', 'S25-2-0007', 'S25-2-0008', 'S25-2-0009', 'S25-2-0010', 'S25-2-0011', 'S25-2-0012', 'S25-2-0013', 'S25-2-0014', 'S25-2-0015', 'S25-2-0016', 'S25-2-0017', 'S25-2-0018', 'S25-2-0019', 'S25-2-0020', 'S25-2-0021', 'S25-2-0022', 'S25-2-0023', 'S25-2-0024', 'S25-2-0025', 'S25-2-0026', 'S25-2-0027', 'S25-2-0028', 'S25-2-0029', 'S25-2-0030', 'S25-2-0031', 'S25-2-0032', 'S25-2-0033', 'S25-2-0034', 'S25-2-0035', 'S25-2-0036', 'S25-2-0037', 'S25-2-0038', 'S25-2-0039', 'S25-2-0040', 'S25-2-0041', 'S25-2-0042', 'S25-2-0043', 'S25-2-0044', 'S25-2-0045', 'S25-2-0046', 'S25-2-0047', 'S25-2-0048', 'S25-2-0049', 'S25-2-0050', 'S25-2-0051', 'S25-2-0052', 'S25-2-0053', 'S25-2-0054', 'S25-2-0055', 'S25-2-0056', 'S25-2-0057', 'S25-2-0058', 'S25-2-0059', 'S25-2-0060', 'S25-2-0061', 'S25-2-0062', 'S25-2-0063'};
G=[ones(1276, 1); -ones(1104, 1)]; % 1 is MDD, -1 controls

% Directories
roisignals='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/Rest-MDD_data/ROISignals_FunImgARCWF/';
motionparams='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/Rest-MDD_data/RealignParameter/';
outputdir='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/matrices/';

stack=nan(58, 58, length(subs));
incomplete=zeros(length(subs), 1);
fd_vols_all=nan(length(subs), 1);

fdthresh=0.25;

% Regions labeled as DMN by Power
dmnidx=[74 75 76 77 78 79 80 81 82 83 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 137 139];

for sub=1:length(subs)
    
    subs{sub}
    
    % load rois data
    load(strcat(roisignals, 'ROISignals_', subs{sub}, '.mat'));
    roifinal=ROISignals;
    
    if size(roifinal, 2)<1833  % some subs have less rois, eliminate them
        disp(strcat(subs{sub}, {' '}, 'data incomplete'))
        incomplete(sub)=1;
    else
        
        % load motion parameters
        fd=dlmread(strcat(motionparams, subs{sub}, '/FD_Power_', subs{sub},'.txt'));
        fd_vols_all(sub)=sum(fd>fdthresh);
        
        % loading Power 264 ROIs
        powerrois=roifinal(:, [1570:1833]);
        
        % Extracting only DMN
        dmnrois=powerrois(:, dmnidx);
        
        
        % if more than 25% volumes have motion, eliminate sub
        if (sum(fd>fdthresh)/size(fd, 1))>0.25
            disp(strcat(subs{sub}, {' '}, 'too much motion'))
            incomplete(sub)=1;
        else
            
            % scrub if FD>0.20
            dmnrois(fd>fdthresh,:)=[];
            
            % check if there are nans in the timeseries
            if sum(sum(isnan(dmnrois)))>0
                disp(strcat(subs{sub}, {' '}, 'data incomplete'))
                incomplete(sub)=1;
            else
                
                % compute correlation and Fisher transform
                corrmat=atanh(corr(dmnrois));
                
                % check if correlation is nan
                if sum(sum(isnan(corrmat)))>0
                    disp(strcat(subs{sub}, {' '}, 'correlation is nan'))
                    incomplete(sub)=1;
                else
                    
                    % save matrix in stack
                    stack(:,:,sub)=corrmat;
                    
                end
            end
        end
    end
end

stack(:,:,logical(incomplete))=[];
subs(logical(incomplete))=[];
G(logical(incomplete))=[];
fd_vols_all(logical(incomplete))=[];

% Save data
save(strcat(outputdir,'connmatstack_dmn_metamdd'), 'stack');
save(strcat(outputdir,'subslist_dmn_metamdd'),'subs')
save(strcat(outputdir,'group_dmn_metamdd'),'G')
save(strcat(outputdir,'fd_volumes_dmn_metamdd'), 'fd_vols_all')





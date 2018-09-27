import pyxiannotator

XiAnnotator = pyxiannotator.XiAnnotator.XiAnnotatorLocal()


peptide = pyxiannotator.AnnotatedSpectrum.Peptide(
    pep_seq1='NIKNDGSDVIVATGFNQQEALQAISSDDRR',
    pep_seq2='YLADKNDLAKVGFIFVDGQIEK',
    link_pos1=3,
    link_pos2=10,
)


peak_list = [[367.609405518, 3292.73], [369.236663818, 3432.24], [389.240356445, 4528.38], [410.240570068, 4475.87], [442.212158203, 3415.13], [465.953277588, 4307.09], [466.244415283, 3698.04], [574.321594238, 60705.4], [575.325073242, 11601.3], [596.596191406, 4780.4], [616.322021484, 4666.06], [689.350280762, 27241.2], [690.353820801, 5523.99], [710.959594727, 4936.0], [730.315856934, 3759.09], [738.67010498, 3772.57], [745.387390137, 5188.45], [788.41607666, 27684.5], [808.906921387, 6828.38], [809.408508301, 7428.29], [820.387634277, 12714.6], [821.394836426, 4627.37], [846.44934082, 5966.95], [846.788452148, 6724.32], [847.117675781, 5932.85], [847.460876465, 4506.28], [857.440307617, 4957.14], [857.773986816, 8708.52], [865.927978516, 27401.3], [866.432495117, 36896.3], [866.931640625, 13076.4], [867.429748535, 6089.5], [919.469116211, 5246.14], [935.488647461, 59409.1], [936.486755371, 30063.4], [939.46105957, 95062.5], [939.964660645, 85812.4], [940.464904785, 46673.4], [940.971130371, 10910.2], [967.973205566, 199214.0], [968.475097656, 229009.0], [968.976257324, 104187.0], [969.477966309, 39791.0], [969.979553223, 10120.1], [1004.51135254, 17007.1], [1018.49639893, 121099.0], [1018.99859619, 112462.0], [1019.4987793, 60801.8], [1054.01513672, 123918.0], [1054.51757813, 130188.0], [1055.01855469, 85277.9], [1103.55126953, 168357.0], [1104.05310059, 181517.0], [1105.88415527, 285440.0], [1106.21728516, 553105.0], [1106.55163574, 529950.0], [1106.88537598, 301703.0], [1107.22070313, 152800.0], [1116.87414551, 261592.0], [1117.20812988, 247174.0], [1177.99645996, 214939.0], [1178.19482422, 281531.0], [1178.39428711, 317084.0], [1178.59448242, 333236.0], [1178.79418945, 253296.0], [1181.79736328, 15242900.0], [1181.99438477, 12525300.0], [1182.19006348, 11736800.0], [1182.38574219, 6855140.0], [1182.5859375, 2934000.0], [1182.78796387, 1023500.0], [1210.12353516, 64706.0], [1233.29345703, 66896.6], [1233.62976074, 73507.7], [1252.66113281, 113099.0], [1269.17260742, 352905.0], [1269.67468262, 472571.0], [1270.17553711, 332863.0], [1270.67578125, 159486.0], [1285.15979004, 420059.0], [1285.66064453, 632321.0], [1286.16101074, 417757.0], [1286.6619873, 242547.0], [1287.16137695, 91472.8], [1294.66467285, 110979.0], [1295.16247559, 81318.4], [1323.67724609, 171987.0], [1324.01293945, 186442.0], [1324.34558105, 167987.0], [1324.67675781, 107118.0], [1339.67041016, 116787.0], [1342.68530273, 87741.2], [1343.01806641, 92254.7], [1343.3515625, 85843.4], [1352.72668457, 25490.8], [1379.43798828, 29016.5], [1379.68884277, 33599.7], [1379.94555664, 33227.4], [1380.18347168, 44382.3], [1392.0411377, 32670.8], [1396.67980957, 39873.0], [1397.18676758, 52426.5], [1397.68725586, 40322.2], [1408.21057129, 31300.8], [1439.9720459, 29807.2], [1440.21850586, 60424.1], [1440.46862793, 61389.1], [1440.71936035, 45959.2], [1440.96875, 60530.5], [1445.19335938, 27604.7], [1445.69091797, 29395.2], [1453.70117188, 83945.9], [1454.20202637, 135061.0], [1454.70568848, 76506.9], [1455.20605469, 33837.7], [1521.74694824, 6699.86], [1527.75927734, 6405.35], [1533.85974121, 6482.97], [1534.84448242, 6928.54], [1536.73913574, 7189.42], [1544.76477051, 39032.8], [1545.25793457, 37933.5], [1545.76220703, 23837.5], [1546.26159668, 14503.3], [1546.75878906, 7730.49], [1560.74353027, 12506.0], [1561.24194336, 20451.1], [1561.74633789, 15381.0], [1565.80957031, 7734.73], [1602.84716797, 6818.66], [1604.85913086, 10987.0], [1635.82055664, 8456.25], [1636.84033203, 6642.72], [1637.84228516, 7402.06], [1658.82141113, 24431.6], [1659.32678223, 18466.4], [1659.83337402, 13232.2], [1660.31860352, 7137.53], [1674.30004883, 12916.9], [1674.80334473, 13469.0], [1675.31677246, 15682.3], [1675.80871582, 11030.6], [1683.30822754, 15443.3], [1683.8112793, 17290.3], [1684.31665039, 14552.5], [1684.8145752, 9728.89], [1749.91760254, 5893.69], [1979.01049805, 6623.39], [1996.00061035, 8972.27]]


annotation_request = XiAnnotator.create_json_annotation_request(
    peak_list=peak_list,
    peptide=peptide,
    precursor_charge=5,
    precursor_mz=1181.190615066879,
    fragment_types=['peptide', 'b', 'y'],
    fragment_tolerance_ppm=10,
    cross_linker="DSSO",
    custom_settings=''
)

XiAnnotator.request_annotation_json(annotation_request)
annotated_spectrum = XiAnnotator.get_annotated_spectrum()
 
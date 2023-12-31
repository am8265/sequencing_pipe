//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.RunMetrics {

using Illumina.InterOp.Run;
using Illumina.InterOp.Metrics;

public class c_csharp_run_metrics {
  public static void copy_focus(base_extraction_metrics metrics, float[] focus_scores, uint channel, uint n) {
    unsafe{ fixed ( float* swig_ptrTo_focus_scores = focus_scores ) {
    {
      c_csharp_run_metricsPINVOKE.copy_focus(base_extraction_metrics.getCPtr(metrics), (global::System.IntPtr)swig_ptrTo_focus_scores, channel, n);
      if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    }
    }}
  }

  public static uint count_q_metric_bins(base_q_metrics metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.count_q_metric_bins__SWIG_1(base_q_metrics.getCPtr(metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint count_q_metric_bins(base_q_by_lane_metrics metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.count_q_metric_bins__SWIG_2(base_q_by_lane_metrics.getCPtr(metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint count_q_metric_bins(base_q_collapsed_metrics arg0) {
    uint ret = c_csharp_run_metricsPINVOKE.count_q_metric_bins__SWIG_3(base_q_collapsed_metrics.getCPtr(arg0));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static void populate_cumulative_distribution(base_q_by_lane_metrics q_metric_set) {
    c_csharp_run_metricsPINVOKE.populate_cumulative_distribution__SWIG_0(base_q_by_lane_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void populate_cumulative_distribution(base_q_metrics q_metric_set) {
    c_csharp_run_metricsPINVOKE.populate_cumulative_distribution__SWIG_1(base_q_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void populate_cumulative_distribution(base_q_collapsed_metrics q_metric_set) {
    c_csharp_run_metricsPINVOKE.populate_cumulative_distribution__SWIG_2(base_q_collapsed_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static bool requires_legacy_bins(uint count) {
    bool ret = c_csharp_run_metricsPINVOKE.requires_legacy_bins(count);
    return ret;
  }

  public static void populate_legacy_q_score_bins(q_score_bin_vector q_score_bins, instrument_type instrument, uint count) {
    c_csharp_run_metricsPINVOKE.populate_legacy_q_score_bins__SWIG_0(q_score_bin_vector.getCPtr(q_score_bins), (int)instrument, count);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void compress_q_metrics(base_q_metrics q_metric_set) {
    c_csharp_run_metricsPINVOKE.compress_q_metrics__SWIG_0(base_q_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void compress_q_metrics(base_q_by_lane_metrics q_metric_set) {
    c_csharp_run_metricsPINVOKE.compress_q_metrics__SWIG_1(base_q_by_lane_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void populate_legacy_q_score_bins(base_q_metrics q_metric_set, q_score_bin_vector q_score_bins, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.populate_legacy_q_score_bins__SWIG_1(base_q_metrics.getCPtr(q_metric_set), q_score_bin_vector.getCPtr(q_score_bins), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void populate_legacy_q_score_bins(base_q_by_lane_metrics q_metric_set, q_score_bin_vector q_score_bins, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.populate_legacy_q_score_bins__SWIG_2(base_q_by_lane_metrics.getCPtr(q_metric_set), q_score_bin_vector.getCPtr(q_score_bins), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static uint count_qvals(base_q_metrics q_metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.count_qvals__SWIG_0(base_q_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint count_qvals(base_q_by_lane_metrics q_metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.count_qvals__SWIG_1(base_q_by_lane_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint count_qvals(base_q_collapsed_metrics arg0) {
    uint ret = c_csharp_run_metricsPINVOKE.count_qvals__SWIG_2(base_q_collapsed_metrics.getCPtr(arg0));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static bool is_compressed(base_q_metrics q_metric_set) {
    bool ret = c_csharp_run_metricsPINVOKE.is_compressed__SWIG_0(base_q_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static bool is_compressed(base_q_collapsed_metrics q_metric_set) {
    bool ret = c_csharp_run_metricsPINVOKE.is_compressed__SWIG_1(base_q_collapsed_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static bool is_compressed(base_q_by_lane_metrics q_metric_set) {
    bool ret = c_csharp_run_metricsPINVOKE.is_compressed__SWIG_2(base_q_by_lane_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint max_qval(base_q_metrics q_metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.max_qval__SWIG_0(base_q_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint max_qval(base_q_by_lane_metrics q_metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.max_qval__SWIG_1(base_q_by_lane_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint max_qval(base_q_collapsed_metrics q_metric_set) {
    uint ret = c_csharp_run_metricsPINVOKE.max_qval__SWIG_2(base_q_collapsed_metrics.getCPtr(q_metric_set));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint index_for_q_value(base_q_metrics q_metric_set, uint qval) {
    uint ret = c_csharp_run_metricsPINVOKE.index_for_q_value__SWIG_0(base_q_metrics.getCPtr(q_metric_set), qval);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static uint index_for_q_value(base_q_by_lane_metrics q_metric_set, uint qval) {
    uint ret = c_csharp_run_metricsPINVOKE.index_for_q_value__SWIG_1(base_q_by_lane_metrics.getCPtr(q_metric_set), qval);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public static void create_collapse_q_metrics(base_q_metrics metric_set, base_q_collapsed_metrics collapsed) {
    c_csharp_run_metricsPINVOKE.create_collapse_q_metrics(base_q_metrics.getCPtr(metric_set), base_q_collapsed_metrics.getCPtr(collapsed));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void create_q_metrics_by_lane(base_q_metrics metric_set, base_q_by_lane_metrics bylane, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.create_q_metrics_by_lane(base_q_metrics.getCPtr(metric_set), base_q_by_lane_metrics.getCPtr(bylane), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static metric_group to_group(metric_type type) {
    metric_group ret = (metric_group)c_csharp_run_metricsPINVOKE.to_group((int)type);
    return ret;
  }

  public static string to_description(metric_type type) {
    string ret = c_csharp_run_metricsPINVOKE.to_description((int)type);
    return ret;
  }

  public static void list_descriptions(metric_type_description_vector types) {
    c_csharp_run_metricsPINVOKE.list_descriptions(metric_type_description_vector.getCPtr(types));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static metric_feature_type to_feature(metric_type type) {
    metric_feature_type ret = (metric_feature_type)c_csharp_run_metricsPINVOKE.to_feature__SWIG_0((int)type);
    return ret;
  }

  public static metric_feature_type to_feature(metric_group type) {
    metric_feature_type ret = (metric_feature_type)c_csharp_run_metricsPINVOKE.to_feature__SWIG_1((int)type);
    return ret;
  }

  public static bool is_base_metric(metric_type type) {
    bool ret = c_csharp_run_metricsPINVOKE.is_base_metric((int)type);
    return ret;
  }

  public static bool is_channel_metric(metric_type type) {
    bool ret = c_csharp_run_metricsPINVOKE.is_channel_metric((int)type);
    return ret;
  }

  public static bool is_read_metric(metric_type type) {
    bool ret = c_csharp_run_metricsPINVOKE.is_read_metric((int)type);
    return ret;
  }

  public static bool is_cycle_metric(metric_type type) {
    bool ret = c_csharp_run_metricsPINVOKE.is_cycle_metric((int)type);
    return ret;
  }

  public static bool is_tile_metric(metric_type type) {
    bool ret = c_csharp_run_metricsPINVOKE.is_tile_metric((int)type);
    return ret;
  }

  public static void list_metrics_to_load(metric_group group, uchar_vector valid_to_load, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_0((int)group, uchar_vector.getCPtr(valid_to_load), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_group group, uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_1((int)group, uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_type type, uchar_vector valid_to_load, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_2((int)type, uchar_vector.getCPtr(valid_to_load), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_type type, uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_3((int)type, uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_group_vector groups, uchar_vector valid_to_load, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_4(metric_group_vector.getCPtr(groups), uchar_vector.getCPtr(valid_to_load), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_group_vector groups, uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_5(metric_group_vector.getCPtr(groups), uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_type_vector types, uchar_vector valid_to_load, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_6(metric_type_vector.getCPtr(types), uchar_vector.getCPtr(valid_to_load), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(metric_type_vector types, uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_7(metric_type_vector.getCPtr(types), uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(string metric_name, uchar_vector valid_to_load, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_8(metric_name, uchar_vector.getCPtr(valid_to_load), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_metrics_to_load(string metric_name, uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_metrics_to_load__SWIG_9(metric_name, uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_summary_metric_groups(metric_group_vector groups, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_summary_metric_groups__SWIG_0(metric_group_vector.getCPtr(groups), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_summary_metric_groups(metric_group_vector groups) {
    c_csharp_run_metricsPINVOKE.list_summary_metric_groups__SWIG_1(metric_group_vector.getCPtr(groups));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_index_summary_metric_groups(metric_group_vector groups) {
    c_csharp_run_metricsPINVOKE.list_index_summary_metric_groups(metric_group_vector.getCPtr(groups));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_summary_metrics_to_load(uchar_vector valid_to_load, instrument_type instrument) {
    c_csharp_run_metricsPINVOKE.list_summary_metrics_to_load__SWIG_0(uchar_vector.getCPtr(valid_to_load), (int)instrument);
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_summary_metrics_to_load(uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_summary_metrics_to_load__SWIG_1(uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_index_metrics_to_load(uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_index_metrics_to_load(uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public static void list_analysis_metrics_to_load(uchar_vector valid_to_load) {
    c_csharp_run_metricsPINVOKE.list_analysis_metrics_to_load(uchar_vector.getCPtr(valid_to_load));
    if (c_csharp_run_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_run_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

}

}

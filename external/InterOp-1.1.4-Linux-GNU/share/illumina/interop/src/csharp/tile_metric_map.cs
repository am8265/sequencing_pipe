//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------

namespace Illumina.InterOp.Metrics {

using System;
using System.Runtime.InteropServices;
using Illumina.InterOp.Run;

public class tile_metric_map : global::System.IDisposable 
    , global::System.Collections.Generic.IDictionary<ulong, base_metric>
 {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal tile_metric_map(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(tile_metric_map obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~tile_metric_map() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          c_csharp_metricsPINVOKE.delete_tile_metric_map(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }


  public base_metric this[ulong key] {
    get {
      return getitem(key);
    }

    set {
      setitem(key, value);
    }
  }

  public bool TryGetValue(ulong key, out base_metric value) {
    if (this.ContainsKey(key)) {
      value = this[key];
      return true;
    }
    value = default(base_metric);
    return false;
  }

  public int Count {
    get {
      return (int)size();
    }
  }

  public bool IsReadOnly {
    get { 
      return false; 
    }
  }

  public global::System.Collections.Generic.ICollection<ulong> Keys {
    get {
      global::System.Collections.Generic.ICollection<ulong> keys = new global::System.Collections.Generic.List<ulong>();
      int size = this.Count;
      if (size > 0) {
        global::System.IntPtr iter = create_iterator_begin();
        for (int i = 0; i < size; i++) {
          keys.Add(get_next_key(iter));
        }
        destroy_iterator(iter);
      }
      return keys;
    }
  }

  public global::System.Collections.Generic.ICollection<base_metric> Values {
    get {
      global::System.Collections.Generic.ICollection<base_metric> vals = new global::System.Collections.Generic.List<base_metric>();
      foreach (global::System.Collections.Generic.KeyValuePair<ulong, base_metric> pair in this) {
        vals.Add(pair.Value);
      }
      return vals;
    }
  }
  
  public void Add(global::System.Collections.Generic.KeyValuePair<ulong, base_metric> item) {
    Add(item.Key, item.Value);
  }

  public bool Remove(global::System.Collections.Generic.KeyValuePair<ulong, base_metric> item) {
    if (Contains(item)) {
      return Remove(item.Key);
    } else {
      return false;
    }
  }

  public bool Contains(global::System.Collections.Generic.KeyValuePair<ulong, base_metric> item) {
    if (this[item.Key] == item.Value) {
      return true;
    } else {
      return false;
    }
  }

  public void CopyTo(global::System.Collections.Generic.KeyValuePair<ulong, base_metric>[] array) {
    CopyTo(array, 0);
  }

  public void CopyTo(global::System.Collections.Generic.KeyValuePair<ulong, base_metric>[] array, int arrayIndex) {
    if (array == null)
      throw new global::System.ArgumentNullException("array");
    if (arrayIndex < 0)
      throw new global::System.ArgumentOutOfRangeException("arrayIndex", "Value is less than zero");
    if (array.Rank > 1)
      throw new global::System.ArgumentException("Multi dimensional array.", "array");
    if (arrayIndex+this.Count > array.Length)
      throw new global::System.ArgumentException("Number of elements to copy is too large.");

    global::System.Collections.Generic.IList<ulong> keyList = new global::System.Collections.Generic.List<ulong>(this.Keys);
    for (int i = 0; i < keyList.Count; i++) {
      ulong currentKey = keyList[i];
      array.SetValue(new global::System.Collections.Generic.KeyValuePair<ulong, base_metric>(currentKey, this[currentKey]), arrayIndex+i);
    }
  }

  global::System.Collections.Generic.IEnumerator<global::System.Collections.Generic.KeyValuePair<ulong, base_metric>> global::System.Collections.Generic.IEnumerable<global::System.Collections.Generic.KeyValuePair<ulong, base_metric>>.GetEnumerator() {
    return new tile_metric_mapEnumerator(this);
  }

  global::System.Collections.IEnumerator global::System.Collections.IEnumerable.GetEnumerator() {
    return new tile_metric_mapEnumerator(this);
  }

  public tile_metric_mapEnumerator GetEnumerator() {
    return new tile_metric_mapEnumerator(this);
  }

  // Type-safe enumerator
  /// Note that the IEnumerator documentation requires an InvalidOperationException to be thrown
  /// whenever the collection is modified. This has been done for changes in the size of the
  /// collection but not when one of the elements of the collection is modified as it is a bit
  /// tricky to detect unmanaged code that modifies the collection under our feet.
  public sealed class tile_metric_mapEnumerator : global::System.Collections.IEnumerator, 
      global::System.Collections.Generic.IEnumerator<global::System.Collections.Generic.KeyValuePair<ulong, base_metric>>
  {
    private tile_metric_map collectionRef;
    private global::System.Collections.Generic.IList<ulong> keyCollection;
    private int currentIndex;
    private object currentObject;
    private int currentSize;

    public tile_metric_mapEnumerator(tile_metric_map collection) {
      collectionRef = collection;
      keyCollection = new global::System.Collections.Generic.List<ulong>(collection.Keys);
      currentIndex = -1;
      currentObject = null;
      currentSize = collectionRef.Count;
    }

    // Type-safe iterator Current
    public global::System.Collections.Generic.KeyValuePair<ulong, base_metric> Current {
      get {
        if (currentIndex == -1)
          throw new global::System.InvalidOperationException("Enumeration not started.");
        if (currentIndex > currentSize - 1)
          throw new global::System.InvalidOperationException("Enumeration finished.");
        if (currentObject == null)
          throw new global::System.InvalidOperationException("Collection modified.");
        return (global::System.Collections.Generic.KeyValuePair<ulong, base_metric>)currentObject;
      }
    }

    // Type-unsafe IEnumerator.Current
    object global::System.Collections.IEnumerator.Current {
      get {
        return Current;
      }
    }

    public bool MoveNext() {
      int size = collectionRef.Count;
      bool moveOkay = (currentIndex+1 < size) && (size == currentSize);
      if (moveOkay) {
        currentIndex++;
        ulong currentKey = keyCollection[currentIndex];
        currentObject = new global::System.Collections.Generic.KeyValuePair<ulong, base_metric>(currentKey, collectionRef[currentKey]);
      } else {
        currentObject = null;
      }
      return moveOkay;
    }

    public void Reset() {
      currentIndex = -1;
      currentObject = null;
      if (collectionRef.Count != currentSize) {
        throw new global::System.InvalidOperationException("Collection modified.");
      }
    }

    public void Dispose() {
      currentIndex = -1;
      currentObject = null;
    }
  }
  

  public tile_metric_map() : this(c_csharp_metricsPINVOKE.new_tile_metric_map__SWIG_0(), true) {
  }

  public tile_metric_map(tile_metric_map other) : this(c_csharp_metricsPINVOKE.new_tile_metric_map__SWIG_1(tile_metric_map.getCPtr(other)), true) {
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  private uint size() {
    uint ret = c_csharp_metricsPINVOKE.tile_metric_map_size(swigCPtr);
    return ret;
  }

  public bool empty() {
    bool ret = c_csharp_metricsPINVOKE.tile_metric_map_empty(swigCPtr);
    return ret;
  }

  public void Clear() {
    c_csharp_metricsPINVOKE.tile_metric_map_Clear(swigCPtr);
  }

  private base_metric getitem(ulong key) {
    base_metric ret = new base_metric(c_csharp_metricsPINVOKE.tile_metric_map_getitem(swigCPtr, key), false);
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private void setitem(ulong key, base_metric x) {
    c_csharp_metricsPINVOKE.tile_metric_map_setitem(swigCPtr, key, base_metric.getCPtr(x));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public bool ContainsKey(ulong key) {
    bool ret = c_csharp_metricsPINVOKE.tile_metric_map_ContainsKey(swigCPtr, key);
    return ret;
  }

  public void Add(ulong key, base_metric val) {
    c_csharp_metricsPINVOKE.tile_metric_map_Add(swigCPtr, key, base_metric.getCPtr(val));
    if (c_csharp_metricsPINVOKE.SWIGPendingException.Pending) throw c_csharp_metricsPINVOKE.SWIGPendingException.Retrieve();
  }

  public bool Remove(ulong key) {
    bool ret = c_csharp_metricsPINVOKE.tile_metric_map_Remove(swigCPtr, key);
    return ret;
  }

  private global::System.IntPtr create_iterator_begin() {
    global::System.IntPtr ret = c_csharp_metricsPINVOKE.tile_metric_map_create_iterator_begin(swigCPtr);
    return ret;
  }

  private ulong get_next_key(global::System.IntPtr swigiterator) {
    ulong ret = c_csharp_metricsPINVOKE.tile_metric_map_get_next_key(swigCPtr, swigiterator);
    return ret;
    }

  private void destroy_iterator(global::System.IntPtr swigiterator) {
    c_csharp_metricsPINVOKE.tile_metric_map_destroy_iterator(swigCPtr, swigiterator);
  }

}

}

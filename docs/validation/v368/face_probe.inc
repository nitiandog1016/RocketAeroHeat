// 仅读取本轮真实组装中指定五条边的面值；活跃通道各写一行。
if constexpr (nDim == 2) if (TailBudget::capturing) {
  for (size_t k = 0; k < Double::Size; ++k) {
    if (!updateMask[k]) continue;
    const auto a = geometry.nodes->GetGlobalIndex(iPoint[k]);
    const auto b = geometry.nodes->GetGlobalIndex(jPoint[k]);
    if (std::min(a,b) < 35*109+10 || std::max(a,b) > 35*109+15 ||
        std::max(a,b)-std::min(a,b) != 1) continue;
    #pragma omp critical(m030_face_probe)
    {
      FILE* f = std::fopen("face_samples.csv", "ab");
      if (!f) std::abort();
      std::fprintf(f, "%lu,%lu,%lu", (unsigned long)iEdge[k], a, b);
      auto put = [&](double x) { std::fprintf(f, ",%.17g", x); };
      put(normal(0)[k]); put(normal(1)[k]);
      put(nonPhysical[k]); put(umusclRamp);
      for (size_t j = 0; j < 6; ++j) put(V1st.i.all(j)[k]);
      for (size_t j = 0; j < 6; ++j) put(V1st.j.all(j)[k]);
      for (size_t j = 0; j < 6; ++j) put(V.i.all(j)[k]);
      for (size_t j = 0; j < 6; ++j) put(V.j.all(j)[k]);
      for (size_t j = 0; j < 4; ++j) put(solution.GetLimiter_Primitive(iPoint[k], j));
      for (size_t j = 0; j < 4; ++j) put(solution.GetLimiter_Primitive(jPoint[k], j));
      put(centralMassSample[k]); put(TailBudget::central[k]);
      put(flux(0)[k]); put(flux(3)[k]);
      std::fputc('\n', f);
      if (std::fclose(f)) std::abort();
    }
  }
}

cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/<file>
declare_datasource (FILE exp_01.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/exp_01.fasta
                    URL_HASH SHA256=e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06)
declare_datasource (FILE IBF_1
                    URL ${CMAKE_SOURCE_DIR}/test/data/IBF_1
                    URL_HASH SHA256=c41b34c5853e2ca738d89171f2ca04662d99e7f0d7498b0f5c555cdc9af7d814)
declare_datasource (FILE mini_example.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.fasta
                    URL_HASH SHA256=f872221632e7423b071d1aedaf4c4e9da2a659e843fcafac12cd65632b904b93)
declare_datasource (FILE mini_example.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.header
                    URL_HASH SHA256=577bcf1c9d018abf60416fac281c4995fa7e300b83cf2d404336fe64507a3bbb)
declare_datasource (FILE mini_example.minimiser
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example.minimiser
                    URL_HASH SHA256=98927a56465368db8c8a557f0ad5b83c25f57fb1820ba997be23b69fb4fe9244)
declare_datasource (FILE mini_example2.header
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example2.header
                    URL_HASH SHA256=9df57d832dea68fe2ef19d783337f9193521b3ad3664b6431d79cd807a362e9f)
declare_datasource (FILE mini_gen.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen.fasta
                    URL_HASH SHA256=6e9da2f6693938586c902f5e4445f2df1f1ac94cff8c23dea9e02b58759a8998)
declare_datasource (FILE mini_gen2.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_gen2.fasta
                    URL_HASH SHA256=7e21b3eff20f950e6a0fa2369a0bbbed4cc967a24083876350f839f4a232771b)
declare_datasource (FILE mini_genom.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_genom.fasta
                    URL_HASH SHA256=b53541bb438be5241880828048fd37544c6aca30d363db46d0d918a2bc531a0e)

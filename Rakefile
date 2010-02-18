
require 'rubygems'
require 'rake'
require 'jeweler'
require 'rake/testtask'
# require 'rcov/rcovtask'

NAME = "dna_sequence_aligner"
WEBSITE_BASE = "website"
WEBSITE_OUTPUT = WEBSITE_BASE + "/output"

gemspec = Gem::Specification.new do |s|
  s.name = NAME
  s.authors = ["John T. Prince"]
  s.email = "jtprince@gmail.com"
  s.homepage = "http://jtprince.github.com/" + NAME
  s.summary = "does high pairwise alignment of sequencing reads with a template"
  s.description = "does high pairwise alignment of sequencing reads with a template using bioruby and clustalw.  gives template-centric output."
  s.add_dependency("bio")
  s.add_development_dependency("spec-more")
end

Jeweler::Tasks.new(gemspec)

Rake::TestTask.new(:spec) do |spec|
  spec.libs << 'lib' << 'spec'
  spec.pattern = 'spec/**/*_spec.rb'
  spec.verbose = true
end

require 'rake/rdoctask'
Rake::RDocTask.new do |rdoc|
  base_rdoc_output_dir = WEBSITE_OUTPUT + '/rdoc'
  version = File.read('VERSION')
  rdoc.rdoc_dir = base_rdoc_output_dir + "/#{version}"
  rdoc.title = NAME + ' ' + version
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end

task :default => :spec

task :build => :gemspec

# credit: Rakefile modeled after Jeweler's

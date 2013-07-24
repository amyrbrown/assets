


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>Redirector/core/model/LPSolver.py at master · bionomicron/Redirector · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://github.global.ssl.fastly.net/images/modules/logos_page/Octocat.png">
    <meta name="hostname" content="fe4.rs.github.com">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="xhr-socket" href="/_sockets" />
    
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" />

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="DaC8d6fe8Rckxy6bshSE97q+Tc/v+yPDi951fkuaWz0=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-ad8bba1e22f69fcea2f2822cf877987d50a76547.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-02d8290450626963b8341bcf949ab6f840ed92b3.css" media="all" rel="stylesheet" type="text/css" />
    


      <script src="https://github.global.ssl.fastly.net/assets/frameworks-e8054ad804a1cf9e9849130fee5a4a5487b663ed.js" type="text/javascript"></script>
      <script src="https://github.global.ssl.fastly.net/assets/github-a9c9533b5e2305f95bfe7370e447c1c68a68897b.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="2713bbdda7640aef08d34d02c218315e">

        <link data-pjax-transient rel='permalink' href='/bionomicron/Redirector/blob/302e8cff531a58d5061977d807fd60838f4e3988/core/model/LPSolver.py'>
  <meta property="og:title" content="Redirector"/>
  <meta property="og:type" content="githubog:gitrepository"/>
  <meta property="og:url" content="https://github.com/bionomicron/Redirector"/>
  <meta property="og:image" content="https://github.global.ssl.fastly.net/images/gravatars/gravatar-user-420.png"/>
  <meta property="og:site_name" content="GitHub"/>
  <meta property="og:description" content="Redirector framework for finding metabolic engineering targets to optimize production and other microbiological phenotypes"/>

  <meta name="description" content="Redirector framework for finding metabolic engineering targets to optimize production and other microbiological phenotypes" />

  <meta content="3010465" name="octolytics-dimension-user_id" /><meta content="bionomicron" name="octolytics-dimension-user_login" /><meta content="7099615" name="octolytics-dimension-repository_id" /><meta content="bionomicron/Redirector" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="7099615" name="octolytics-dimension-repository_network_root_id" /><meta content="bionomicron/Redirector" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/bionomicron/Redirector/commits/master.atom" rel="alternate" title="Recent Commits to Redirector:master" type="application/atom+xml" />

  </head>


  <body class="logged_out page-blob  vis-public env-production ">

    <div class="wrapper">
      
      
      


      
      <div class="header header-logged-out">
  <div class="container clearfix">

    <a class="header-logo-wordmark" href="https://github.com/">
      <span class="mega-octicon octicon-logo-github"></span>
    </a>

    <div class="header-actions">
        <a class="button primary" href="/signup">Sign up</a>
      <a class="button" href="/login?return_to=%2Fbionomicron%2FRedirector%2Fblob%2Fmaster%2Fcore%2Fmodel%2FLPSolver.py">Sign in</a>
    </div>

    <div class="command-bar js-command-bar  in-repository">


      <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
        <li class="features"><a href="/features">Features</a></li>
          <li class="enterprise"><a href="https://enterprise.github.com/">Enterprise</a></li>
          <li class="blog"><a href="/blog">Blog</a></li>
      </ul>
        <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    
      data-repo="bionomicron/Redirector"
      data-branch="master"
      data-sha="54a44e1c672406b92dd3b0b61f86b09d73f80a2e"
  >

    <input type="hidden" name="nwo" value="bionomicron/Redirector" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
    </div>

  </div>
</div>


      


          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">


  <li>
  <a href="/login?return_to=%2Fbionomicron%2FRedirector"
  class="minibutton with-count js-toggler-target star-button entice tooltipped upwards "
  title="You must be signed in to use this feature" rel="nofollow">
  <span class="octicon octicon-star"></span>Star
</a>
<a class="social-count js-social-count" href="/bionomicron/Redirector/stargazers">
  0
</a>

  </li>

    <li>
      <a href="/login?return_to=%2Fbionomicron%2FRedirector"
        class="minibutton with-count js-toggler-target fork-button entice tooltipped upwards"
        title="You must be signed in to fork a repository" rel="nofollow">
        <span class="octicon octicon-git-branch"></span>Fork
      </a>
      <a href="/bionomicron/Redirector/network" class="social-count">
        0
      </a>
    </li>
</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author">
            <a href="/bionomicron" class="url fn" itemprop="url" rel="author"><span itemprop="title">bionomicron</span></a></span
          ><span class="repohead-name-divider">/</span><strong
          ><a href="/bionomicron/Redirector" class="js-current-repository js-repo-home-link">Redirector</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">

      <div class="repository-with-sidebar repo-container ">

        <div class="repository-sidebar">
            

<div class="repo-nav repo-nav-full js-repository-container-pjax js-octicon-loaders">
  <div class="repo-nav-contents">
    <ul class="repo-menu">
      <li class="tooltipped leftwards" title="Code">
        <a href="/bionomicron/Redirector" aria-label="Code" class="js-selected-navigation-item selected" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /bionomicron/Redirector">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped leftwards" title="Issues">
          <a href="/bionomicron/Redirector/issues" aria-label="Issues" class="js-selected-navigation-item js-disable-pjax" data-gotokey="i" data-selected-links="repo_issues /bionomicron/Redirector/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped leftwards" title="Pull Requests"><a href="/bionomicron/Redirector/pulls" aria-label="Pull Requests" class="js-selected-navigation-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /bionomicron/Redirector/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/bionomicron/Redirector/wiki" aria-label="Wiki" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_wiki /bionomicron/Redirector/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="repo-menu-separator"></div>
    <ul class="repo-menu">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/bionomicron/Redirector/pulse" aria-label="Pulse" class="js-selected-navigation-item " data-pjax="true" data-selected-links="pulse /bionomicron/Redirector/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/bionomicron/Redirector/graphs" aria-label="Graphs" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_graphs repo_contributors /bionomicron/Redirector/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/bionomicron/Redirector/network" aria-label="Network" class="js-selected-navigation-item js-disable-pjax" data-selected-links="repo_network /bionomicron/Redirector/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

    </ul>

  </div>
</div>

            <div class="only-with-full-nav">
              

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=clone">
  <h3><strong>HTTPS</strong> clone URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/bionomicron/Redirector.git" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/bionomicron/Redirector.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=clone">
  <h3><strong>Subversion</strong> checkout URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/bionomicron/Redirector" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/bionomicron/Redirector" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>



<p class="clone-options">You can clone with
    <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
    <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>,
  and <a href="https://help.github.com/articles/which-remote-url-should-i-use">other methods.</a>
</p>



                <a href="/bionomicron/Redirector/archive/master.zip"
                   class="minibutton sidebar-button"
                   title="Download this repository as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
            </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:c9dd8f6392ca12dcb5a3fdd38004449f -->
<!-- blob contrib frag key: views10/v8/blob_contributors:v21:c9dd8f6392ca12dcb5a3fdd38004449f -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/bionomicron/Redirector/find/master" data-pjax data-hotkey="t" style="display:none">Show File Finder</a>

<div class="file-navigation">
  


<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/bionomicron/Redirector/blob/master/core/model/LPSolver.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="master" data-skip-pjax="true" rel="nofollow" title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/bionomicron/Redirector" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">Redirector</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/bionomicron/Redirector/tree/master/core" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">core</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/bionomicron/Redirector/tree/master/core/model" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">model</span></a></span><span class="separator"> / </span><strong class="final-path">LPSolver.py</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="core/model/LPSolver.py" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


  
  <div class="commit file-history-tease">
    <img class="main-avatar" height="24" src="https://secure.gravatar.com/avatar/1f083079606ab10ba7573a0f74aed70b?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png" width="24" />
    <span class="author"><a href="/bionomicron" rel="author">bionomicron</a></span>
    <time class="js-relative-date" datetime="2013-04-02T11:24:00-07:00" title="2013-04-02 11:24:00">April 02, 2013</time>
    <div class="commit-title">
        <a href="/bionomicron/Redirector/commit/c63be814c9b46de79acefe1d05b9058b8d782a8f" class="message" data-pjax="true" title="fixed issue with parsing gene names">fixed issue with parsing gene names</a>
    </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong> contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
        <li class="facebox-user-list-item">
          <img height="24" src="https://secure.gravatar.com/avatar/1f083079606ab10ba7573a0f74aed70b?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png" width="24" />
          <a href="/bionomicron">bionomicron</a>
        </li>
      </ul>
    </div>
  </div>


<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
          <span>334 lines (267 sloc)</span>
        <span>10.765 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
              <a class="minibutton js-entice" href=""
                 data-entice="You must be signed in and on a branch to make or propose changes">Edit</a>
          <a href="/bionomicron/Redirector/raw/master/core/model/LPSolver.py" class="button minibutton " id="raw-url">Raw</a>
            <a href="/bionomicron/Redirector/blame/master/core/model/LPSolver.py" class="button minibutton ">Blame</a>
          <a href="/bionomicron/Redirector/commits/master/core/model/LPSolver.py" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
            <a class="minibutton danger empty-icon js-entice" href=""
               data-entice="You must be signed in and on a branch to make or propose changes">
            Delete
          </a>
      </div><!-- /.actions -->

    </div>
        <div class="blob-wrapper data type-python js-blob-data">
      <table class="file-code file-diff">
        <tr class="file-code-line">
          <td class="blob-line-nums">
            <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>
<span id="L287" rel="#L287">287</span>
<span id="L288" rel="#L288">288</span>
<span id="L289" rel="#L289">289</span>
<span id="L290" rel="#L290">290</span>
<span id="L291" rel="#L291">291</span>
<span id="L292" rel="#L292">292</span>
<span id="L293" rel="#L293">293</span>
<span id="L294" rel="#L294">294</span>
<span id="L295" rel="#L295">295</span>
<span id="L296" rel="#L296">296</span>
<span id="L297" rel="#L297">297</span>
<span id="L298" rel="#L298">298</span>
<span id="L299" rel="#L299">299</span>
<span id="L300" rel="#L300">300</span>
<span id="L301" rel="#L301">301</span>
<span id="L302" rel="#L302">302</span>
<span id="L303" rel="#L303">303</span>
<span id="L304" rel="#L304">304</span>
<span id="L305" rel="#L305">305</span>
<span id="L306" rel="#L306">306</span>
<span id="L307" rel="#L307">307</span>
<span id="L308" rel="#L308">308</span>
<span id="L309" rel="#L309">309</span>
<span id="L310" rel="#L310">310</span>
<span id="L311" rel="#L311">311</span>
<span id="L312" rel="#L312">312</span>
<span id="L313" rel="#L313">313</span>
<span id="L314" rel="#L314">314</span>
<span id="L315" rel="#L315">315</span>
<span id="L316" rel="#L316">316</span>
<span id="L317" rel="#L317">317</span>
<span id="L318" rel="#L318">318</span>
<span id="L319" rel="#L319">319</span>
<span id="L320" rel="#L320">320</span>
<span id="L321" rel="#L321">321</span>
<span id="L322" rel="#L322">322</span>
<span id="L323" rel="#L323">323</span>
<span id="L324" rel="#L324">324</span>
<span id="L325" rel="#L325">325</span>
<span id="L326" rel="#L326">326</span>
<span id="L327" rel="#L327">327</span>
<span id="L328" rel="#L328">328</span>
<span id="L329" rel="#L329">329</span>
<span id="L330" rel="#L330">330</span>
<span id="L331" rel="#L331">331</span>
<span id="L332" rel="#L332">332</span>
<span id="L333" rel="#L333">333</span>
<span id="L334" rel="#L334">334</span>

          </td>
          <td class="blob-line-code">
                  <div class="highlight"><pre><div class='line' id='LC1'><span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC2'><span class="sd">@author: Graham Rockwell</span></div><div class='line' id='LC3'><span class="sd">Updated 2012.05.28</span></div><div class='line' id='LC4'><span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC5'><br/></div><div class='line' id='LC6'><span class="kn">from</span> <span class="nn">pulp</span> <span class="kn">import</span> <span class="n">LpVariable</span><span class="p">,</span> <span class="n">LpAffineExpression</span><span class="p">,</span> <span class="n">LpConstraint</span><span class="p">,</span> <span class="n">LpProblem</span><span class="p">,</span> <span class="n">LpMinimize</span><span class="p">,</span> <span class="n">LpStatus</span></div><div class='line' id='LC7'><span class="kn">import</span> <span class="nn">pulp</span></div><div class='line' id='LC8'><span class="kn">import</span> <span class="nn">os</span><span class="o">,</span><span class="nn">re</span></div><div class='line' id='LC9'><span class="kn">import</span> <span class="nn">tempfile</span></div><div class='line' id='LC10'><br/></div><div class='line' id='LC11'><br/></div><div class='line' id='LC12'><span class="c">#=================================================</span></div><div class='line' id='LC13'><span class="c">#=============LinearOptimization==================</span></div><div class='line' id='LC14'><span class="c">#=================================================</span></div><div class='line' id='LC15'><br/></div><div class='line' id='LC16'><span class="k">class</span> <span class="nc">LPSolver</span><span class="p">:</span></div><div class='line' id='LC17'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC18'><span class="sd">    LP solver class</span></div><div class='line' id='LC19'><span class="sd">    Used for solving LP and MILP problems</span></div><div class='line' id='LC20'><span class="sd">    &#39;&#39;&#39;</span></div><div class='line' id='LC21'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC22'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">__init__</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">name</span><span class="o">=</span><span class="s">&#39;analysis&#39;</span><span class="p">,</span><span class="n">verbose</span><span class="o">=</span><span class="bp">False</span><span class="p">):</span></div><div class='line' id='LC23'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span> <span class="o">=</span> <span class="n">LpProblem</span><span class="p">()</span></div><div class='line' id='LC24'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpVariables</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC25'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpObjective</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC26'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC27'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">objectiveValue</span> <span class="o">=</span> <span class="bp">None</span></div><div class='line' id='LC28'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">rowNames</span><span class="o">=</span> <span class="p">[]</span></div><div class='line' id='LC30'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">columnNames</span><span class="o">=</span> <span class="p">[]</span></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">configFile</span> <span class="o">=</span> <span class="s">&#39;&#39;</span></div><div class='line' id='LC32'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">statusCode</span> <span class="o">=</span> <span class="bp">None</span></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">verbose</span> <span class="o">=</span> <span class="n">verbose</span></div><div class='line' id='LC35'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">useLimitTag</span> <span class="o">=</span> <span class="bp">False</span></div><div class='line' id='LC36'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">mpsLogFile</span> <span class="o">=</span> <span class="bp">False</span></div><div class='line' id='LC37'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">ignoreBadReferences</span> <span class="o">=</span> <span class="bp">False</span></div><div class='line' id='LC38'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">Mip</span> <span class="o">=</span> <span class="bp">True</span></div><div class='line' id='LC39'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC40'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">_parseSCIPLog</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">fileName</span><span class="p">):</span></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC42'><span class="sd">        Parse result of SCIP analysis.</span></div><div class='line' id='LC43'><span class="sd">        This section needs a much more solid testing as it key interface with the solver </span></div><div class='line' id='LC44'><span class="sd">        but often solver errors are dropped when they occur with no / minimal warnings.</span></div><div class='line' id='LC45'><span class="sd">        </span></div><div class='line' id='LC46'><span class="sd">        @var fileName: name and directory information of file being parsed, usually a temp file </span></div><div class='line' id='LC47'><span class="sd">        @type fileName: string</span></div><div class='line' id='LC48'><span class="sd">        @return: (values of fluxes, value of objective function, status of optimization)</span></div><div class='line' id='LC49'><span class="sd">        @rtype: (dict,float,string) </span></div><div class='line' id='LC50'><span class="sd">        &#39;&#39;&#39;</span></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fileHandle</span> <span class="o">=</span> <span class="nb">open</span><span class="p">(</span><span class="n">fileName</span><span class="p">,</span> <span class="s">&#39;r&#39;</span><span class="p">)</span></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">location</span> <span class="o">=</span> <span class="mi">0</span></div><div class='line' id='LC53'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC54'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC55'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">while</span> <span class="ow">not</span> <span class="n">line</span><span class="o">.</span><span class="n">startswith</span><span class="p">(</span><span class="s">&#39;SCIP Status        :&#39;</span><span class="p">):</span></div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">line</span><span class="o">.</span><span class="n">startswith</span><span class="p">(</span><span class="s">&#39;Syntax error&#39;</span><span class="p">):</span></div><div class='line' id='LC57'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span> <span class="s">&quot;SCIP parsing syntax error: </span><span class="si">%s</span><span class="s">&quot;</span> <span class="o">%</span> <span class="n">line</span></div><div class='line' id='LC58'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC59'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span> <span class="n">line</span></div><div class='line' id='LC60'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ilocation</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">tell</span><span class="p">()</span> </div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span>  <span class="n">location</span> <span class="o">==</span> <span class="n">ilocation</span><span class="p">:</span></div><div class='line' id='LC62'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span> <span class="s">&quot;failed to read file </span><span class="si">%s</span><span class="s">&quot;</span> <span class="o">%</span><span class="p">(</span><span class="n">fileName</span><span class="p">)</span></div><div class='line' id='LC63'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="p">({},</span><span class="mf">0.0</span><span class="p">,</span><span class="s">&quot;failed to read file&quot;</span><span class="p">)</span></div><div class='line' id='LC64'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">location</span> <span class="o">=</span> <span class="n">ilocation</span></div><div class='line' id='LC65'><br/></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">m</span> <span class="o">=</span> <span class="n">re</span><span class="o">.</span><span class="n">match</span><span class="p">(</span><span class="s">&#39;[\w\s:]*\[([\w\s]*)\]&#39;</span><span class="p">,</span> <span class="n">line</span><span class="p">)</span></div><div class='line' id='LC69'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">statusString</span> <span class="o">=</span> <span class="n">m</span><span class="o">.</span><span class="n">group</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span></div><div class='line' id='LC70'><br/></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">objectiveValue</span> <span class="o">=</span> <span class="bp">None</span></div><div class='line' id='LC72'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">=</span> <span class="nb">dict</span><span class="p">()</span></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">statusString</span> <span class="o">==</span> <span class="s">&#39;optimal solution found&#39;</span><span class="p">:</span></div><div class='line' id='LC75'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC76'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">while</span> <span class="ow">not</span> <span class="n">line</span><span class="o">.</span><span class="n">startswith</span><span class="p">(</span><span class="s">&#39;objective value:&#39;</span><span class="p">):</span></div><div class='line' id='LC77'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c">#Good place to check all the features of the solution.</span></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">line</span><span class="o">.</span><span class="n">startswith</span><span class="p">(</span><span class="s">&#39;solution violates&#39;</span><span class="p">):</span></div><div class='line' id='LC79'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span> <span class="s">&quot;Error solving optimization [</span><span class="si">%s</span><span class="s">]&quot;</span> <span class="o">%</span> <span class="n">line</span></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ilocation</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">tell</span><span class="p">()</span> </div><div class='line' id='LC81'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span>  <span class="n">location</span> <span class="o">==</span> <span class="n">ilocation</span><span class="p">:</span></div><div class='line' id='LC82'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span> <span class="s">&quot;failed to read file </span><span class="si">%s</span><span class="s">&quot;</span> <span class="o">%</span><span class="p">(</span><span class="n">fileName</span><span class="p">)</span></div><div class='line' id='LC83'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="p">({},</span><span class="mf">0.0</span><span class="p">,</span><span class="s">&quot;failed to read file&quot;</span><span class="p">)</span></div><div class='line' id='LC84'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC85'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC86'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">m</span> <span class="o">=</span> <span class="n">re</span><span class="o">.</span><span class="n">match</span><span class="p">(</span><span class="s">&#39;objective value:\s+([\S]+)&#39;</span><span class="p">,</span> <span class="n">line</span><span class="p">)</span></div><div class='line' id='LC87'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">objectiveValue</span> <span class="o">=</span> <span class="nb">float</span><span class="p">(</span><span class="n">m</span><span class="o">.</span><span class="n">group</span><span class="p">(</span><span class="mi">1</span><span class="p">))</span></div><div class='line' id='LC88'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC90'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">while</span> <span class="ow">not</span> <span class="n">line</span><span class="o">.</span><span class="n">strip</span><span class="p">()</span> <span class="o">==</span> <span class="s">&#39;&#39;</span><span class="p">:</span></div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">m</span> <span class="o">=</span> <span class="n">re</span><span class="o">.</span><span class="n">match</span><span class="p">(</span><span class="s">&#39;([\S]+)\s+([\S]+)&#39;</span><span class="p">,</span> <span class="n">line</span><span class="p">)</span></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">name</span> <span class="o">=</span> <span class="n">m</span><span class="o">.</span><span class="n">group</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span></div><div class='line' id='LC93'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span><span class="p">[</span><span class="n">name</span><span class="p">]</span> <span class="o">=</span> <span class="nb">float</span><span class="p">(</span><span class="n">m</span><span class="o">.</span><span class="n">group</span><span class="p">(</span><span class="mi">2</span><span class="p">))</span></div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">line</span> <span class="o">=</span> <span class="n">fileHandle</span><span class="o">.</span><span class="n">readline</span><span class="p">()</span></div><div class='line' id='LC95'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fileHandle</span><span class="o">.</span><span class="n">close</span><span class="p">()</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC98'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="p">(</span><span class="n">result</span><span class="p">,</span> <span class="n">objectiveValue</span><span class="p">,</span> <span class="n">statusString</span><span class="p">)</span></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC100'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">writeMPS</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">fileName</span><span class="p">):</span></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&quot;&quot;&quot;</span></div><div class='line' id='LC103'><span class="sd">        @param fileName: name of file</span></div><div class='line' id='LC104'><span class="sd">        @type fileName: String</span></div><div class='line' id='LC105'><span class="sd">        &quot;&quot;&quot;</span></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span><span class="o">.</span><span class="n">writeMPS</span><span class="p">(</span><span class="n">fileName</span><span class="p">)</span></div><div class='line' id='LC107'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">True</span></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC109'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">writeLP</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">fileName</span><span class="p">):</span></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&quot;&quot;&quot;</span></div><div class='line' id='LC111'><span class="sd">        @param fileName: name of file</span></div><div class='line' id='LC112'><span class="sd">        @type fileName: String</span></div><div class='line' id='LC113'><span class="sd">        &quot;&quot;&quot;</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span><span class="o">.</span><span class="n">writeLP</span><span class="p">(</span><span class="n">fileName</span><span class="p">,</span><span class="n">writeSOS</span><span class="o">=</span><span class="mi">0</span><span class="p">)</span></div><div class='line' id='LC115'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">True</span></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC117'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">clear</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC119'><span class="sd">        Left over place holder function</span></div><div class='line' id='LC120'><span class="sd">        To be removed</span></div><div class='line' id='LC121'><span class="sd">        &#39;&#39;&#39;</span></div><div class='line' id='LC122'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">None</span></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC124'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">setConfigFile</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">configFile</span><span class="p">):</span></div><div class='line' id='LC125'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC126'><span class="sd">        @param configFile: name of configuration file handed to solver (mostly for SCIP)</span></div><div class='line' id='LC127'><span class="sd">        @type configFile: String</span></div><div class='line' id='LC128'><span class="sd">        &#39;&#39;&#39;</span></div><div class='line' id='LC129'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">configFile</span> <span class="o">=</span> <span class="n">configFile</span></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">setObjectiveByName</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">name</span><span class="p">,</span><span class="n">value</span><span class="p">):</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC133'><span class="sd">        @param name: name of variable to be optimized</span></div><div class='line' id='LC134'><span class="sd">        @param value: coefficient of objective </span></div><div class='line' id='LC135'><span class="sd">        </span></div><div class='line' id='LC136'><span class="sd">        Useful for setting a optimization of a single variable as the objective function </span></div><div class='line' id='LC137'><span class="sd">        &#39;&#39;&#39;</span></div><div class='line' id='LC138'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpObjective</span> <span class="o">=</span> <span class="p">{</span><span class="n">name</span><span class="p">:</span><span class="n">value</span><span class="p">}</span></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC140'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">setObjectiveFunction</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">objectiveMap</span><span class="p">):</span></div><div class='line' id='LC141'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpObjective</span> <span class="o">=</span> <span class="n">objectiveMap</span></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC143'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">clearObjective</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpObjective</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC145'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">_getLpVariables</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">model</span><span class="p">):</span></div><div class='line' id='LC147'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">columnNames</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getColumnNames</span><span class="p">()</span></div><div class='line' id='LC148'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">columnLimits</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getColumnLimits</span><span class="p">()</span></div><div class='line' id='LC149'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">mipColumns</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getMipColumnNames</span><span class="p">()</span></div><div class='line' id='LC150'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC151'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC152'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">varNames</span> <span class="o">=</span> <span class="p">[]</span></div><div class='line' id='LC153'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC154'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="n">name</span> <span class="ow">in</span> <span class="n">columnNames</span><span class="p">:</span></div><div class='line' id='LC155'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">(</span><span class="n">lower</span><span class="p">,</span><span class="n">upper</span><span class="p">)</span> <span class="o">=</span> <span class="p">(</span><span class="bp">None</span><span class="p">,</span><span class="bp">None</span><span class="p">)</span></div><div class='line' id='LC156'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">name</span> <span class="ow">in</span> <span class="n">columnLimits</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></div><div class='line' id='LC157'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">(</span><span class="n">lower</span><span class="p">,</span><span class="n">upper</span><span class="p">)</span> <span class="o">=</span> <span class="n">columnLimits</span><span class="p">[</span><span class="n">name</span><span class="p">]</span></div><div class='line' id='LC158'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="nb">type</span> <span class="o">=</span> <span class="s">&quot;Continuous&quot;</span></div><div class='line' id='LC159'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">name</span> <span class="ow">in</span> <span class="n">mipColumns</span><span class="p">:</span></div><div class='line' id='LC160'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="nb">type</span> <span class="o">=</span> <span class="s">&quot;Integer&quot;</span></div><div class='line' id='LC161'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span><span class="p">[</span><span class="n">name</span><span class="p">]</span> <span class="o">=</span> <span class="n">LpVariable</span><span class="p">(</span><span class="n">name</span><span class="p">,</span><span class="n">lower</span><span class="p">,</span><span class="n">upper</span><span class="p">,</span><span class="n">cat</span><span class="o">=</span><span class="nb">type</span><span class="p">)</span></div><div class='line' id='LC162'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC163'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">result</span></div><div class='line' id='LC164'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC165'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">_getLpProblem</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">model</span><span class="p">,</span><span class="n">lpVariables</span><span class="p">,</span><span class="n">objective</span><span class="o">=</span><span class="bp">None</span><span class="p">):</span></div><div class='line' id='LC166'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">modelName</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">modelName</span></div><div class='line' id='LC167'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">rowLimits</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getRowLimits</span><span class="p">()</span></div><div class='line' id='LC168'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">objective</span> <span class="o">==</span> <span class="bp">None</span><span class="p">:</span></div><div class='line' id='LC169'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">objective</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getObjective</span><span class="p">()</span></div><div class='line' id='LC170'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC171'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">exp</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC172'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">=</span> <span class="n">LpProblem</span><span class="p">(</span><span class="n">modelName</span><span class="p">,</span><span class="n">LpMinimize</span><span class="p">)</span></div><div class='line' id='LC173'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC174'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">lpObjective</span> <span class="o">=</span> <span class="n">LpAffineExpression</span><span class="p">()</span></div><div class='line' id='LC175'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="n">key</span><span class="p">,</span><span class="n">value</span><span class="p">)</span> <span class="ow">in</span> <span class="n">objective</span><span class="o">.</span><span class="n">items</span><span class="p">():</span></div><div class='line' id='LC176'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">lpObjective</span> <span class="o">+=</span> <span class="n">lpVariables</span><span class="p">[</span><span class="n">key</span><span class="p">]</span><span class="o">*</span><span class="n">value</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC178'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpObjective</span> <span class="o">=</span> <span class="n">lpObjective</span></div><div class='line' id='LC179'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">+=</span> <span class="n">lpObjective</span></div><div class='line' id='LC180'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC181'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="n">key</span><span class="p">,</span><span class="n">value</span><span class="p">)</span> <span class="ow">in</span> <span class="n">model</span><span class="o">.</span><span class="n">data</span><span class="o">.</span><span class="n">data</span><span class="o">.</span><span class="n">items</span><span class="p">():</span></div><div class='line' id='LC182'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">(</span><span class="n">rowName</span><span class="p">,</span><span class="n">columnName</span><span class="p">)</span> <span class="o">=</span> <span class="n">key</span></div><div class='line' id='LC183'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC184'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">rowName</span> <span class="ow">not</span> <span class="ow">in</span> <span class="n">exp</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></div><div class='line' id='LC185'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">exp</span><span class="p">[</span><span class="n">rowName</span><span class="p">]</span> <span class="o">=</span> <span class="n">LpAffineExpression</span><span class="p">()</span></div><div class='line' id='LC186'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">exp</span><span class="p">[</span><span class="n">rowName</span><span class="p">]</span> <span class="o">+=</span> <span class="n">lpVariables</span><span class="p">[</span><span class="n">columnName</span><span class="p">]</span><span class="o">*</span><span class="n">value</span></div><div class='line' id='LC187'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC188'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="n">key</span><span class="p">,(</span><span class="n">lower</span><span class="p">,</span><span class="n">upper</span><span class="p">))</span> <span class="ow">in</span> <span class="n">rowLimits</span><span class="o">.</span><span class="n">items</span><span class="p">():</span></div><div class='line' id='LC189'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">key</span> <span class="ow">not</span> <span class="ow">in</span> <span class="n">exp</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></div><div class='line' id='LC190'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">continue</span></div><div class='line' id='LC191'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC192'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">tupper</span> <span class="o">=</span> <span class="n">upper</span> <span class="ow">not</span> <span class="ow">in</span> <span class="p">[</span><span class="bp">None</span><span class="p">,</span><span class="nb">float</span><span class="p">(</span><span class="s">&quot;inf&quot;</span><span class="p">)]</span></div><div class='line' id='LC193'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">tlower</span> <span class="o">=</span> <span class="n">lower</span> <span class="ow">not</span> <span class="ow">in</span> <span class="p">[</span><span class="bp">None</span><span class="p">,</span><span class="nb">float</span><span class="p">(</span><span class="s">&quot;-inf&quot;</span><span class="p">)]</span></div><div class='line' id='LC194'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">tag</span> <span class="o">=</span> <span class="s">&#39;&#39;</span></div><div class='line' id='LC195'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC196'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">lower</span> <span class="o">==</span> <span class="n">upper</span> <span class="ow">and</span> <span class="n">tupper</span> <span class="ow">and</span> <span class="n">tlower</span><span class="p">:</span></div><div class='line' id='LC197'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">useLimitTag</span><span class="p">:</span> <span class="n">tag</span> <span class="o">=</span> <span class="s">&quot;</span><span class="si">%s</span><span class="s">&quot;</span> <span class="o">%</span> <span class="p">(</span><span class="n">key</span><span class="p">)</span></div><div class='line' id='LC198'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">+=</span> <span class="n">exp</span><span class="p">[</span><span class="n">key</span><span class="p">]</span> <span class="o">==</span> <span class="n">upper</span><span class="p">,</span> <span class="n">tag</span></div><div class='line' id='LC199'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC200'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">tupper</span><span class="p">:</span></div><div class='line' id='LC201'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">useLimitTag</span><span class="p">:</span> <span class="n">tag</span> <span class="o">=</span> <span class="s">&quot;</span><span class="si">%s</span><span class="s">_u&quot;</span> <span class="o">%</span> <span class="p">(</span><span class="n">key</span><span class="p">)</span></div><div class='line' id='LC202'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">+=</span> <span class="n">exp</span><span class="p">[</span><span class="n">key</span><span class="p">]</span> <span class="o">&lt;=</span> <span class="n">upper</span><span class="p">,</span> <span class="n">tag</span></div><div class='line' id='LC203'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">tlower</span><span class="p">:</span></div><div class='line' id='LC204'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">useLimitTag</span><span class="p">:</span> <span class="n">tag</span> <span class="o">=</span> <span class="s">&quot;</span><span class="si">%s</span><span class="s">_l&quot;</span> <span class="o">%</span> <span class="p">(</span><span class="n">key</span><span class="p">)</span></div><div class='line' id='LC205'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">+=</span> <span class="n">exp</span><span class="p">[</span><span class="n">key</span><span class="p">]</span> <span class="o">&gt;=</span> <span class="n">lower</span><span class="p">,</span> <span class="n">tag</span></div><div class='line' id='LC206'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="ow">not</span> <span class="n">tupper</span> <span class="ow">and</span> <span class="ow">not</span> <span class="n">tlower</span><span class="p">:</span></div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">pass</span></div><div class='line' id='LC208'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC209'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">result</span></div><div class='line' id='LC210'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC211'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">setModel</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">model</span><span class="p">,</span><span class="n">objective</span> <span class="o">=</span> <span class="bp">None</span><span class="p">):</span></div><div class='line' id='LC212'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">columnNames</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getColumnNames</span><span class="p">()</span></div><div class='line' id='LC213'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">rowNames</span> <span class="o">=</span> <span class="n">model</span><span class="o">.</span><span class="n">getRowNames</span><span class="p">()</span></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC215'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpVariables</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">_getLpVariables</span><span class="p">(</span><span class="n">model</span><span class="p">)</span></div><div class='line' id='LC216'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">_getLpProblem</span><span class="p">(</span><span class="n">model</span><span class="p">,</span> <span class="bp">self</span><span class="o">.</span><span class="n">lpVariables</span><span class="p">,</span> <span class="n">objective</span><span class="p">)</span></div><div class='line' id='LC217'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC218'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="p">(</span><span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span><span class="p">,</span><span class="bp">self</span><span class="o">.</span><span class="n">lpVariables</span><span class="p">)</span></div><div class='line' id='LC219'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC220'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">_testLP</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">testKeys</span><span class="p">):</span></div><div class='line' id='LC221'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">c</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span><span class="o">.</span><span class="n">constraints</span></div><div class='line' id='LC222'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ks</span> <span class="o">=</span> <span class="n">c</span><span class="o">.</span><span class="n">keys</span><span class="p">()</span></div><div class='line' id='LC223'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ks</span><span class="o">.</span><span class="n">sort</span><span class="p">()</span></div><div class='line' id='LC224'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="n">k</span> <span class="ow">in</span> <span class="n">ks</span><span class="p">:</span></div><div class='line' id='LC225'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">constraint</span> <span class="o">=</span> <span class="n">c</span><span class="p">[</span><span class="n">k</span><span class="p">]</span></div><div class='line' id='LC226'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="ow">not</span> <span class="n">constraint</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></div><div class='line' id='LC227'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c">#empty constraint add the dummyVar</span></div><div class='line' id='LC228'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">constraint</span> <span class="o">+=</span> <span class="bp">self</span><span class="o">.</span><span class="n">get_dummyVar</span><span class="p">()</span></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ctag</span> <span class="o">=</span> <span class="n">constraint</span><span class="o">.</span><span class="n">asCplexLpConstraint</span><span class="p">(</span><span class="n">k</span><span class="p">)</span></div><div class='line' id='LC230'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">k</span> <span class="ow">in</span> <span class="n">testKeys</span><span class="p">:</span></div><div class='line' id='LC231'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span> <span class="s">&quot;[</span><span class="si">%s</span><span class="s">]: </span><span class="si">%s</span><span class="s">&quot;</span> <span class="o">%</span> <span class="p">(</span><span class="n">k</span><span class="p">,</span><span class="n">ctag</span><span class="p">)</span></div><div class='line' id='LC232'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC233'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">None</span></div><div class='line' id='LC234'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC235'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">_setPredictionMap</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC237'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="n">key</span><span class="p">,</span><span class="n">value</span><span class="p">)</span> <span class="ow">in</span> <span class="bp">self</span><span class="o">.</span><span class="n">lpVariables</span><span class="o">.</span><span class="n">items</span><span class="p">():</span></div><div class='line' id='LC238'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span><span class="p">[</span><span class="n">key</span><span class="p">]</span> <span class="o">=</span> <span class="n">pulp</span><span class="o">.</span><span class="n">value</span><span class="p">(</span><span class="n">value</span><span class="p">)</span></div><div class='line' id='LC239'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC240'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span> <span class="o">=</span> <span class="n">result</span></div><div class='line' id='LC241'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC242'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">solve</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">solver</span><span class="o">=</span><span class="n">pulp</span><span class="o">.</span><span class="n">GLPK</span><span class="p">(</span><span class="n">msg</span><span class="o">=</span><span class="mi">0</span><span class="p">)):</span></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">status</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span><span class="o">.</span><span class="n">solve</span><span class="p">(</span><span class="n">solver</span><span class="p">)</span></div><div class='line' id='LC244'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">_setPredictionMap</span><span class="p">()</span>    </div><div class='line' id='LC245'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">status</span></div><div class='line' id='LC246'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC247'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">runGLPK</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span><span class="n">solver</span><span class="o">=</span><span class="n">pulp</span><span class="o">.</span><span class="n">GLPK</span><span class="p">(</span><span class="n">msg</span><span class="o">=</span><span class="mi">0</span><span class="p">)):</span></div><div class='line' id='LC248'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">status</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">lpProblem</span><span class="o">.</span><span class="n">solve</span><span class="p">(</span><span class="n">solver</span><span class="p">)</span></div><div class='line' id='LC249'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">_setPredictionMap</span><span class="p">()</span>    </div><div class='line' id='LC250'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">status</span></div><div class='line' id='LC251'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC252'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">runSCIP</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC253'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&quot;&quot;&quot;</span></div><div class='line' id='LC254'><span class="sd">        Performs Simplex or MILP solving</span></div><div class='line' id='LC255'><span class="sd">        Calls SCIP solver </span></div><div class='line' id='LC256'><span class="sd">        &quot;&quot;&quot;</span></div><div class='line' id='LC257'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC258'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">Mip</span><span class="o">=</span><span class="bp">True</span></div><div class='line' id='LC259'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">lpTempFile</span> <span class="o">=</span> <span class="n">tempfile</span><span class="o">.</span><span class="n">mkstemp</span><span class="p">(</span><span class="s">&#39;.lp&#39;</span><span class="p">)</span></div><div class='line' id='LC260'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scipTempFile</span> <span class="o">=</span> <span class="n">tempfile</span><span class="o">.</span><span class="n">mkstemp</span><span class="p">()</span></div><div class='line' id='LC261'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">os</span><span class="o">.</span><span class="n">close</span><span class="p">(</span><span class="n">lpTempFile</span><span class="p">[</span><span class="mi">0</span><span class="p">])</span></div><div class='line' id='LC262'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">os</span><span class="o">.</span><span class="n">close</span><span class="p">(</span><span class="n">scipTempFile</span><span class="p">[</span><span class="mi">0</span><span class="p">])</span></div><div class='line' id='LC263'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC264'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">writeLP</span><span class="p">(</span><span class="n">lpTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span></div><div class='line' id='LC265'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC266'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">config</span> <span class="o">=</span> <span class="s">&#39;&#39;</span></div><div class='line' id='LC267'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">configFile</span> <span class="ow">and</span> <span class="bp">self</span><span class="o">.</span><span class="n">configFile</span> <span class="o">!=</span> <span class="s">&#39;&#39;</span><span class="p">:</span></div><div class='line' id='LC268'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c">#config = &quot;-s &quot; + self.configFile</span></div><div class='line' id='LC269'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">pass</span></div><div class='line' id='LC270'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">command</span> <span class="o">=</span> <span class="s">&#39;scip -q -f &#39;</span> <span class="o">+</span> <span class="n">lpTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span> <span class="o">+</span> <span class="s">&#39; -l &#39;</span> <span class="o">+</span> <span class="n">scipTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span> <span class="o">+</span> <span class="s">&quot; &quot;</span><span class="o">+</span> <span class="n">config</span></div><div class='line' id='LC271'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">verbose</span><span class="p">:</span> <span class="k">print</span> <span class="s">&quot;wrote temp file </span><span class="si">%s</span><span class="s">&quot;</span> <span class="o">%</span> <span class="p">(</span><span class="n">scipTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span></div><div class='line' id='LC272'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC273'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">os</span><span class="o">.</span><span class="n">system</span><span class="p">(</span><span class="n">command</span><span class="p">)</span></div><div class='line' id='LC274'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">os</span><span class="o">.</span><span class="n">remove</span><span class="p">(</span><span class="n">lpTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span></div><div class='line' id='LC275'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC276'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">(</span><span class="n">result</span><span class="p">,</span> <span class="bp">self</span><span class="o">.</span><span class="n">objectiveValue</span><span class="p">,</span> <span class="n">status</span><span class="p">)</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">_parseSCIPLog</span><span class="p">(</span><span class="n">scipTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span></div><div class='line' id='LC277'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c">#Need to add check for solution that violates boundaries.            </span></div><div class='line' id='LC278'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">os</span><span class="o">.</span><span class="n">remove</span><span class="p">(</span><span class="n">scipTempFile</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span></div><div class='line' id='LC279'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC280'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span> <span class="o">=</span> <span class="p">{}</span></div><div class='line' id='LC281'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="n">name</span> <span class="ow">in</span> <span class="bp">self</span><span class="o">.</span><span class="n">columnNames</span><span class="p">:</span></div><div class='line' id='LC282'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dName</span> <span class="o">=</span> <span class="n">name</span><span class="o">.</span><span class="n">replace</span><span class="p">(</span><span class="s">&quot;[&quot;</span><span class="p">,</span><span class="s">&quot;(&quot;</span><span class="p">)</span></div><div class='line' id='LC283'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dName</span> <span class="o">=</span> <span class="n">dName</span><span class="o">.</span><span class="n">replace</span><span class="p">(</span><span class="s">&quot;]&quot;</span><span class="p">,</span><span class="s">&quot;)&quot;</span><span class="p">)</span></div><div class='line' id='LC284'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dName</span> <span class="o">=</span> <span class="n">dName</span><span class="o">.</span><span class="n">replace</span><span class="p">(</span><span class="s">&quot;-&quot;</span><span class="p">,</span><span class="s">&quot;~&quot;</span><span class="p">)</span></div><div class='line' id='LC285'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">name</span> <span class="ow">in</span> <span class="n">result</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></div><div class='line' id='LC286'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span><span class="p">[</span><span class="n">name</span><span class="p">]</span> <span class="o">=</span> <span class="n">result</span><span class="p">[</span><span class="n">name</span><span class="p">]</span></div><div class='line' id='LC287'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">elif</span> <span class="n">dName</span> <span class="ow">in</span> <span class="n">result</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></div><div class='line' id='LC288'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span><span class="p">[</span><span class="n">name</span><span class="p">]</span> <span class="o">=</span> <span class="n">result</span><span class="p">[</span><span class="n">dName</span><span class="p">]</span></div><div class='line' id='LC289'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC290'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span><span class="p">[</span><span class="n">name</span><span class="p">]</span> <span class="o">=</span> <span class="mf">0.0</span></div><div class='line' id='LC291'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC292'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c">#status = self.statusMap[reportCode]</span></div><div class='line' id='LC293'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC294'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">status</span></div><div class='line' id='LC295'><br/></div><div class='line' id='LC296'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">runIntOpt</span><span class="p">(</span><span class="bp">self</span><span class="p">,):</span></div><div class='line' id='LC297'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">Mip</span><span class="p">:</span></div><div class='line' id='LC298'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">resultCode</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">runSCIP</span><span class="p">()</span></div><div class='line' id='LC299'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC300'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">resultCode</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">solve</span><span class="p">()</span></div><div class='line' id='LC301'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">resultCode</span></div><div class='line' id='LC302'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC303'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">run</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">model</span> <span class="o">=</span> <span class="bp">None</span><span class="p">,</span> <span class="n">objective</span> <span class="o">=</span> <span class="bp">None</span><span class="p">,</span> <span class="n">logFile</span><span class="o">=</span><span class="bp">None</span><span class="p">):</span></div><div class='line' id='LC304'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="sd">&#39;&#39;&#39;</span></div><div class='line' id='LC305'><span class="sd">        choose and run appropriate optimization</span></div><div class='line' id='LC306'><span class="sd">        &#39;&#39;&#39;</span></div><div class='line' id='LC307'><br/></div><div class='line' id='LC308'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">model</span> <span class="o">!=</span> <span class="bp">None</span><span class="p">:</span></div><div class='line' id='LC309'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">setModel</span><span class="p">(</span><span class="n">model</span><span class="p">,</span><span class="n">objective</span><span class="p">)</span></div><div class='line' id='LC310'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">logFile</span> <span class="o">!=</span> <span class="bp">None</span><span class="p">:</span></div><div class='line' id='LC311'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">writeLP</span><span class="p">(</span><span class="n">logFile</span><span class="p">)</span> </div><div class='line' id='LC312'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="ow">not</span> <span class="bp">self</span><span class="o">.</span><span class="n">Mip</span><span class="p">:</span></div><div class='line' id='LC313'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">resultCode</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">solve</span><span class="p">()</span></div><div class='line' id='LC314'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC315'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">resultCode</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">runSCIP</span><span class="p">()</span></div><div class='line' id='LC316'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">=</span> <span class="bp">self</span><span class="o">.</span><span class="n">getPredictionMap</span><span class="p">()</span></div><div class='line' id='LC317'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC318'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">result</span></div><div class='line' id='LC319'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC320'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">_getLpObjective</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC321'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">result</span> <span class="o">=</span> <span class="n">pulp</span><span class="o">.</span><span class="n">value</span><span class="p">(</span><span class="bp">self</span><span class="o">.</span><span class="n">lpObjective</span><span class="p">)</span></div><div class='line' id='LC322'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">result</span></div><div class='line' id='LC323'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC324'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">getPredictionMap</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC325'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span></div><div class='line' id='LC326'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC327'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">getMipPredictionMap</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC328'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">self</span><span class="o">.</span><span class="n">predictionMap</span></div><div class='line' id='LC329'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC330'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">def</span> <span class="nf">getMipObjectiveValue</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></div><div class='line' id='LC331'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="bp">self</span><span class="o">.</span><span class="n">Mip</span><span class="p">:</span></div><div class='line' id='LC332'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="bp">self</span><span class="o">.</span><span class="n">objectiveValue</span></div><div class='line' id='LC333'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC334'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="bp">self</span><span class="o">.</span><span class="n">_getLpObjective</span><span class="p">()</span></div></pre></div>
          </td>
        </tr>
      </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2013 <span title="0.55055s from fe4.rs.github.com">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/bionomicron/Redirector/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

    
  </body>
</html>


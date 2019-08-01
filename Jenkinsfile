podTemplate(
    name: 'jenkins-portcullis',
    label: 'jenkins-portcullis',
    namespace: 'cicd',
    containers: [
      containerTemplate(name: 'cppbuild', image: 'maplesond/cppbuild:latest', ttyEnabled: true),
      containerTemplate(name: 'gitversion', image: 'maplesond/gitversion:latest', ttyEnabled: true, command: 'cat'),
      containerTemplate(name: 'sonarscanner', image: 'docker.sdlmapleson.net/sonarscanner:1', ttyEnabled: true),
      containerTemplate(name: 'githubrelease', image: 'maplesond/githubrelease:latest', ttyEnabled: true, command: 'cat'),
      containerTemplate(name: 'docker', image: 'docker', ttyEnabled: true, command: 'cat'),
      containerTemplate(name: 'jnlp', image: 'jenkins/jnlp-slave:3.27-1')
    ],
    volumes: [hostPathVolume(mountPath: '/var/run/docker.sock', hostPath: '/var/run/docker.sock')],
    imagePullSecrets: ['docker.sdlmapleson.net'],
    cloud: 'kubernetes') {
  node('jenkins-portcullis') {
    stage('Git') {
      checkout(scm)
    }
    // Only do semantic versioning on master or develop branch for now
    if(env.BRANCH_NAME == 'master' || env.BRANCH_NAME == 'develop') {
      stage('Git Version') {
        container('gitversion') {
          sh ("""SEMVER=`semver .` && echo "\$SEMVER" > version && ./update_version.sh \$SEMVER""")
        }
      }
    }
    // We only have the free sonarqube, so just run on master branch to save time.
    if(env.BRANCH_NAME == 'master') {
      stage('SonarQube') {
        container('sonarscanner') {
          withCredentials([string(credentialsId: 'sonarqube-token', variable: 'TOKEN')]) {
            sh "sonar-scanner -Dsonar.projectKey=portcullis -Dsonar.sources=src,lib,tests,scripts -Dsonar.exclusions=tests/gtest -Dsonar.login=${TOKEN}"
          }
        }
      }
    }    
    stage('Build') {
      container('cppbuild') {
        sh "./autogen.sh"
        sh "./configure"
        sh "make -j4 V=1"
      }
    }
    stage('Test') {
      container('cppbuild') {
        sh "make -j4 V=1 check"
      }
    }
    stage('Install and Package') {
      container('cppbuild') {
        sh "make install"        
        sh "make dist"
      }
    }
    if(env.BRANCH_NAME == 'master' || env.BRANCH_NAME == 'develop') {
      stage('Docker') {
        container('docker') {
          // Get the version we previously saved to disk
          def versionFile = readFile "version"
          SEMVER = versionFile.split('\n')[0]
          // Build the image
          def image = docker.build("harbor.sdlmapleson.net/portcullis/portcullis:${SEMVER}", "--build-arg VERSION=${SEMVER} .")
          image.inside {
            sh "portcullis --help || true"
            sh "junctools --help || true"
          }
          // Push to local registry
          docker.withRegistry('https://harbor.sdlmapleson.net', 'harbor') {
            image.push()
            image.push("latest")
            if(env.BRANCH_NAME == 'master') {
              image.push("stable")
              // Dockerhub
              withCredentials([usernamePassword(credentialsId: 'dockerhub', usernameVariable: 'DOCKERHUB_USER', passwordVariable: 'DOCKERHUB_PASS')]) {
                sh "docker login -u $DOCKERHUB_USER -p $DOCKERHUB_PASS"
                sh "docker tag harbor.sdlmapleson.net/portcullis/portcullis:${SEMVER} maplesond/portcullis:${SEMVER}" 
                sh "docker tag harbor.sdlmapleson.net/portcullis/portcullis:${SEMVER} maplesond/portcullis:latest" 
                sh "docker tag harbor.sdlmapleson.net/portcullis/portcullis:${SEMVER} maplesond/portcullis:stable" 
                sh "docker push maplesond/portcullis:${SEMVER}"
                sh "docker push maplesond/portcullis:latest"
                sh "docker push maplesond/portcullis:stable"
              }
            }
          }
        }
      }
    }
    if(env.BRANCH_NAME == 'master') {
      stage('Github Release') {
        container('githubrelease') {
          withCredentials([string(credentialsId: 'github-maplesond-portcullis', variable: 'GITHUB_TOKEN')]) {
            GITHUB_USER = "maplesond"
            GITHUB_REPO = "portcullis"
            def versionFile = readFile "version"
            SEMVER = versionFile.split('\n')[0]
            sh "git tag ${SEMVER}"
            sh "git remote add github-maplesond https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/maplesond/portcullis.git && git push github-maplesond master"
            sh "git remote add github-ei https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/EI-CoreBioinformatics/portcullis.git && git push github-ei master && git push github-ei master --tags"
            sh """DESCRIPTION=`git log -1 | tail -n +4` && echo "\$DESCRIPTION" > description"""
            def DESCRIPTION = readFile "description"
            sh "ls"
            sh "github-release release -u ${GITHUB_USER} -s ${GITHUB_TOKEN} --tag ${SEMVER} --repo ${GITHUB_REPO} --description '${DESCRIPTION}'"
            sh "github-release upload -u ${GITHUB_USER} -s ${GITHUB_TOKEN} --tag ${SEMVER} --repo ${GITHUB_REPO} --name ${GITHUB_REPO}-${SEMVER}.tar.gz --file ${GITHUB_REPO}-${SEMVER}.tar.gz --label ${GITHUB_REPO}-${SEMVER}.tar.gz"
          }
        }
      }
      // This is probably going to be hard but it should be possible to at least start the process in an automatic way.
      //stage('Brew update') {
      //}
      //stage('Conda update') {
      //}
    }
  }
}
